#pragma once
#include "kinova_dynamics/Model.hpp"
#include "kinova_system/BaseBlock.hpp"
#include "Messages.hpp"
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <deque>

namespace Roahm
{
namespace Dynamics
{
using namespace Eigen;
class ControllerBlock : public ::Roahm::System::BaseBlock
{
  protected:
    using ControlMsg = ::Roahm::msgs::ControlMsg;
    using Measurement = ::Roahm::msgs::Measurement;
    using AccTrajectory = ::Roahm::msgs::AccTrajectory;
    struct InputPack;

  public:
    /**
     * @brief calculation structure of controller
     **/
    virtual void run() final override;

    /**
     * @brief specific function to calculate dynamics
     * @param inputs: %InputPack to use
     **/
    virtual Eigen::VectorXd dynamics(const InputPack &inputs) = 0;

  protected:
    /**
     * @brief Ctor to be called by derived class
     * @param block_name: name of the block
     * @param model_path: path of the model file
     **/
    ControllerBlock(const std::string &block_name,
                    const std::string &model_path);

    /**
     * @brief constructor directly given robot model
     **/
    ControllerBlock(const std::string &block_name,
                    const Model::SpatialModel<double> &rob_model);

    /**
     * @brief package for easy passing of all the value needed
     **/
    struct InputPack
    {
        // desired state
        const VectorXd &pos_des;
        const VectorXd &vel_des;
        const VectorXd &acc_des;

        // actual state
        const VectorXd &pos;
        const VectorXd &vel;

        // pos and vel error
        const VectorXd &e;
        const VectorXd &ed;

        // time difference
        const double delta_t;
    };

  private:
    /**
     * @brief get desired state
     * @param pdes: desired pos
     * @param vdes: desired vel
     * @param ades: desired acc
     * @param t: current timestamp
     **/
    std::shared_ptr<const AccTrajectory> get_desired(
        VectorXd &pdes, VectorXd &vdes, VectorXd &ades,
        decltype(std::chrono::system_clock::now()) t);

  protected:
    /// trajectory acceleration input
    System::AbstractInputPort::SharedPtr traj_port;
    /// state from robot
    System::AbstractInputPort::SharedPtr state_port;

    /// robot model
    Model::SpatialModel<double> rob_model;

    /// time when controller starts working
    decltype(std::chrono::system_clock::now()) t0;

    /// time of the previous message
    decltype(std::chrono::system_clock::now()) t_prev;

    // trajectory queue
    std::deque<std::shared_ptr<const AccTrajectory>> trajectories;
    bool has_trajectory = false;
};
} // namespace Dynamics
} // namespace Roahm
