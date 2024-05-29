#pragma once
#include "kinova_dynamics/ControllerBlock.hpp"
#include "kinova_dynamics/Model.hpp"
#include "kinova_dynamics/RNEA.hpp"

namespace Roahm
{
namespace Dynamics
{
/**
 * @brief RNEA Controller block using PD control
 * @details TODO
 **/
class RNEABlock final : public ControllerBlock
{
  public:
    using SharedPtr = std::shared_ptr<RNEABlock>;

  private:
    using Base = ::Roahm::System::BaseBlock;

  public:
    static SharedPtr make_shared(const std::string &block_name,
                                 const std::string &model_path);

    static SharedPtr make_shared(const std::string &block_name,
                                 const Model::SpatialModel<double> &rob_model);

    /**
     * @brief call RNEA for calculation
     * @param inputs: %InputPack to use
     **/
    virtual Eigen::VectorXd dynamics(const InputPack &inputs) final override;

    /**
     * @brief set coefficients
     * @param K_p
     * @param K_d
     * @param K_r
     **/
    void setK(const MatrixXd &K_p, const MatrixXd &K_d, const MatrixXd &K_r);

  private:
    RNEABlock(const std::string &block_name, const std::string &model_path);
    RNEABlock(const std::string &block_name,
              const Model::SpatialModel<double> &rob_model);

  private:
    // coefficients
    bool is_coeff_set{false};
    MatrixXd K_p;
    MatrixXd K_d;
    MatrixXd K_r;
};
} // namespace Dynamics
} // namespace Roahm
