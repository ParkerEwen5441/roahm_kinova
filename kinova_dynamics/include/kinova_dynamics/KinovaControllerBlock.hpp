#pragma once
#include "kinova_dynamics/ControllerBlock.hpp"
#include "kinova_dynamics/Model.hpp"
#include "kinova_dynamics/RNEA.hpp"

namespace Roahm
{
namespace Dynamics
{
class KinovaControlBlock final : public ControllerBlock
{
  public:
    using SharedPtr = std::shared_ptr<KinovaControlBlock>;

  private:
    using Base = ::Roahm::System::BaseBlock;

  public:
    /**
     * @brief block maker
     **/
    /**
     * @brief Ctor taking non-interval model directly
     * @param block_name: name of control block
     * @param rob_model: model of robot
     * @param variance: variance of the interval
     **/
    KinovaControlBlock(const std::string &block_name,
                       const Model::SpatialModel<double> &rob_model,
                       const double variance, const double fric_var);
    
    static SharedPtr make_shared(const std::string &block_name,
                                 const std::string &model_path,
                                 const double variance,
                                 const double fric_var = -1);

    static SharedPtr make_shared(const std::string &block_name,
                                 const Model::SpatialModel<double> &rob_model,
                                 const double variance,
                                 const double fric_var = -1);

    /**
     * @brief implement robust control
     **/
    virtual Eigen::VectorXd dynamics(const InputPack &inputs) final override;

    /**
     * @brief set coefficients
     * @param K_r
     * @param alpha
     * @param V_max
     * @param r_norm_threshold
     **/
    void setK(const MatrixXd &K_r, const double alpha, const double V_max,
              const double r_norm_threshold);

  private:
    bool is_coeff_set{false};
    double accum_error{0}; // accumulated error
    Model::SpatialModelInt int_model;
    // ARMOUR params
    MatrixXd K_r;
    double alpha;
    double V_max;
    double r_norm_threshold;
};
} // namespace Dynamics
} // namespace Roahm