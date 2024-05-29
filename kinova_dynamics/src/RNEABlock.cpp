#include "kinova_dynamics/RNEABlock.hpp"
#include "kinova_dynamics/ControllerBlock.hpp"
#include "kinova_dynamics/Parser.hpp"
#include "kinova_dynamics/Messages.hpp"
#include <kinova_system/BaseBlock.hpp>
#include <kinova_system/Port/OutputPort.hpp>
#include <stdexcept>

namespace Roahm
{
namespace Dynamics
{

RNEABlock::SharedPtr RNEABlock::make_shared(const std::string &block_name,
                                            const std::string &model_path)
{
    return SharedPtr(new RNEABlock(block_name, model_path));
}

RNEABlock::SharedPtr RNEABlock::make_shared(
    const std::string &block_name, const Model::SpatialModel<double> &rob_model)
{
    return SharedPtr(new RNEABlock(block_name, rob_model));
}

RNEABlock::RNEABlock(const std::string &block_name,
                     const std::string &model_path)
    : ControllerBlock(block_name, model_path)
{
}

RNEABlock::RNEABlock(const std::string &block_name,
                     const Model::SpatialModel<double> &rob_model)
    : ControllerBlock(block_name, rob_model)
{
}

void RNEABlock::setK(const MatrixXd &K_p, const MatrixXd &K_d,
                     const MatrixXd &K_r)
{
    assert(K_p.rows() == rob_model.num_dof);
    assert(K_p.cols() == rob_model.num_dof);
    assert(K_p.rows() == rob_model.num_dof);
    assert(K_d.cols() == rob_model.num_dof);
    assert(K_r.rows() == rob_model.num_dof);
    assert(K_r.cols() == rob_model.num_dof);
    this->K_p = K_p;
    this->K_d = K_d;
    this->K_r = K_r;
    is_coeff_set = true;
}

Eigen::VectorXd RNEABlock::dynamics(const InputPack &inputs)
{
    if (!is_coeff_set)
    {
        const char *err_msg = "No Coefficient for RNEA is set!";
        logger.critical(err_msg);
        throw std::runtime_error(err_msg);
    }
    // VectorXd qd_r = inputs.vel_des + K_r * inputs.e;
    VectorXd qdd_r = inputs.acc_des + K_p * inputs.e + K_d * inputs.ed;
    return rnea(rob_model, inputs.pos, inputs.vel, inputs.vel, qdd_r);
}
} // namespace Dynamics
} // namespace Roahm
