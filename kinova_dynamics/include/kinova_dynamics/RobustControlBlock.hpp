#pragma once
#include "kinova_dynamics/ControllerBlock.hpp"
#include "kinova_dynamics/Model.hpp"
#include "kinova_dynamics/RNEA.hpp"
#include "kinova_system/BaseBlock.hpp"
#include <memory>

namespace Roahm
{
namespace Dynamics
{
class RobustControlBlock final : public ControllerBlock
{
  public:
    using SharedPtr = std::shared_ptr<RobustControlBlock>;

  private:
    using Base = ::Roahm::System::BaseBlock;

  public:
    /**
     * @brief block maker
     **/
    static SharedPtr make_shared(const std::string &block_name,
                                 const std::string &model_path,
                                 const double variance);

    static SharedPtr make_shared(const std::string &block_name,
                                 const Model::SpatialModel<double> &rob_model,
                                 const double variance);

    /**
     * @brief implement robust control
     **/
    virtual Eigen::VectorXd dynamics(const InputPack &inputs) final override;

    /**
     * @brief set coefficients
     * @param K_p
     * @param K_i
     * @param K_r
     * @param maxErrorBound: maximum error bound
     **/
    void setK(const double kappa_p, const double kappa_i, const double phi_p,
              const double phi_i, const MatrixXd &Kr, const double max_error);

  private:
    /**
     * @brief Ctor taking non-interval model directly
     * @param block_name: name of control block
     * @param rob_model: model of robot
     * @param variance: variance of the interval
     **/
    RobustControlBlock(const std::string &block_name,
                       const Model::SpatialModel<double> &rob_model,
                       const double variance);

    bool is_coeff_set{false};
    double accum_error{0}; // accumulated error
    Model::SpatialModelInt int_model;
    MatrixXd K_r;
    double kappa_p;
    double kappa_i;
    double phi_p;
    double phi_i;
    double max_error_bound; // if error > bound, start accumulating
};
} // namespace Dynamics
} // namespace Roahm
