#include "kinova_dynamics/RobustControlBlock.hpp"
#include "kinova_dynamics/Model.hpp"
#include <eigen3/Eigen/src/Core/Matrix.h>

namespace Roahm
{
namespace Dynamics
{
RobustControlBlock::SharedPtr RobustControlBlock::make_shared(
    const std::string &block_name, const std::string &model_path,
    const double variance)
{
    auto rob_model = Model::SpatialModeld::load(model_path);
    return SharedPtr(new RobustControlBlock(block_name, rob_model, variance));
}

RobustControlBlock::SharedPtr RobustControlBlock::make_shared(
    const std::string &block_name, const Model::SpatialModel<double> &rob_model,
    const double variance)
{
    return SharedPtr(new RobustControlBlock(block_name, rob_model, variance));
}

RobustControlBlock::RobustControlBlock(
    const std::string &block_name, const Model::SpatialModel<double> &rob_model,
    const double variance)
    : ControllerBlock(block_name, rob_model),
      int_model(Model::SpatialModelInt::from_double(rob_model, variance))
{
}

Eigen::VectorXd RobustControlBlock::dynamics(const InputPack &inputs)
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
    // error Only do state error for actuated joints
    double state_error = sqrt(inputs.e.dot(inputs.e));
    if (state_error > max_error_bound)
    {
        accum_error += state_error * inputs.delta_t;
    }

    double phi_t = phi_p + phi_i * accum_error;
    double kappa_t = kappa_p + kappa_i * accum_error;

    // Step 7, compute robust input
    VectorXd v = -(kappa_t * bound.norm() + phi_t) * r;

    // Step 8, compute torque
    return u_nominal - v;
}

void RobustControlBlock::setK(const double kappa_p, const double kappa_i,
                              const double phi_p, const double phi_i,
                              const MatrixXd &Kr, const double max_error)
{
    this->kappa_p = kappa_p;
    this->kappa_i = kappa_i;
    this->phi_p = phi_p;
    this->phi_i = phi_i;
    this->max_error_bound = max_error;

    assert(Kr.rows() == rob_model.num_dof);
    assert(Kr.cols() == rob_model.num_dof);
    this->K_r = Kr;

    is_coeff_set = true;
}

} // namespace Dynamics
} // namespace Roahm
