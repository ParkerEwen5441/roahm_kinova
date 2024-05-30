#include "kinova_dynamics/KinovaControllerBlock.hpp"
#include "kinova_dynamics/Model.hpp"
#include "kinova_dynamics/RNEA.hpp"
#include <eigen3/Eigen/src/Core/Matrix.h>

namespace Roahm
{
namespace Dynamics
{
KinovaControlBlock::SharedPtr KinovaControlBlock::make_shared(
    const std::string &block_name, const std::string &model_path,
    const double variance, const double fric_var)
{
    auto rob_model = Model::SpatialModeld::load(model_path);
    return SharedPtr(
        new KinovaControlBlock(block_name, rob_model, variance, fric_var));
}

KinovaControlBlock::SharedPtr KinovaControlBlock::make_shared(
    const std::string &block_name, const Model::SpatialModel<double> &rob_model,
    const double variance, const double fric_var)
{
    return SharedPtr(
        new KinovaControlBlock(block_name, rob_model, variance, fric_var));
}

KinovaControlBlock::KinovaControlBlock(
    const std::string &block_name, const Model::SpatialModel<double> &rob_model,
    const double variance, const double fric_var)
    : ControllerBlock(block_name, rob_model),
      int_model(
          Model::SpatialModelInt::from_double(rob_model, variance, fric_var))
{
}

Eigen::VectorXd KinovaControlBlock::dynamics(const InputPack &inputs)
{
    if (!is_coeff_set)
    {
        const char *err_msg = "Coefficient for Robust Controller isn't set!";
        logger.critical(err_msg);
        throw std::runtime_error(err_msg);
    }

    // Step 1, calculate reference terms
    VectorXd qd_aux = inputs.vel_des + K_r * inputs.e;
    VectorXd qdd_aux = inputs.acc_des + K_r * inputs.ed;
    VectorXd r = inputs.ed + K_r * inputs.e;

    // Step 2, calculate nominal torque from passivity RNEA
    VectorXd u_nominal =
        rnea(rob_model, inputs.pos, inputs.vel, qd_aux, qdd_aux);

    // Step 3, calculate interval torque from passivity RNEA
    VecX<Intd> u_interval =
        int_rnea(int_model, inputs.pos, inputs.vel, qd_aux, qdd_aux);

    // Interval check
    for (size_t i = 0; i < rob_model.num_dof; i++)
        if (u_nominal[i] > u_interval[i].upper ||
            u_nominal[i] < u_interval[i].lower)
            logger.error(
                "Nominal model output falls outside interval output!!!");

    // Step 4, calculate error between nominal and interval
    // Because of constraints, unactuated terms should be 0
    VecX<Intd> phi = u_interval - u_nominal.cast<Intd>();

    // Step 5, calculate error bound (only includes actuated terms)
    VectorXd bound(rob_model.num_dof);
    for (uint32_t i = 0; i < rob_model.num_dof; i++)
    {
        bound[i] = fmax(abs(phi[i].lower), abs(phi[i].upper));
    }

    // Step 6, calculate y(t), k(t), only includes actuated terms for state
    // error

    // Only compute robust input when r is large enough to avoid numerical
    // issues
    VectorXd v = VectorXd::Zero(rob_model.num_dof);
    double r_norm = r.norm();
    if (r_norm > r_norm_threshold)
    {
        // Step 6a, calculate V = 0.5 * r' * M * r
        static const VectorXd zero_input = VectorXd::Zero(rob_model.num_dof);
        VecX<Intd> M_r =
            int_rnea(int_model, inputs.pos, zero_input, zero_input, r, false);

        // dot product between r and M * r
        Intd V_int = 0;
        for (uint32_t i = 0; i < rob_model.num_dof; i++)
        {
            V_int += 0.5 * r(i) * M_r(i);
        }

        // Step 6b, compute robust input
        double h = -V_int.upper + V_max;

        double gamma = fmax(0.0, -alpha * h / r_norm + bound.norm());

        v = gamma * r / r_norm;
    }

    // Step 7, compute torque
    return u_nominal + v;
} // namespace Dynamics

void KinovaControlBlock::setK(const MatrixXd &K_r, const double alpha,
                              const double V_max, const double r_norm_threshold)
{
    this->alpha = alpha;
    this->V_max = V_max;
    this->r_norm_threshold = r_norm_threshold;
    assert(K_r.rows() == rob_model.num_dof);
    assert(K_r.cols() == rob_model.num_dof);
    this->K_r = K_r;

    is_coeff_set = true;
}

} // namespace Dynamics
} // namespace Roahm