#include "kinova_dynamics/ControllerBlock.hpp"
#include "kinova_dynamics/Helper.hpp"
#include "kinova_dynamics/Model.hpp"
#include "kinova_dynamics/Parser.hpp"
#include "kinova_system/BaseBlock.hpp"
#include "kinova_system/Port/InputPort.hpp"
#include "kinova_system/Port/OutputPort.hpp"
#include <chrono>
#include <cmath>
#include <eigen3/Eigen/Eigen>
#include <stdexcept>
using namespace std;

namespace Roahm
{
namespace Dynamics
{
void ControllerBlock::run()
{
    // init output
    ControlMsg ctrl(rob_model.num_dof);
    default_output->set(ctrl);

    // initialization
    // init values
    {
        logger.info("Check input initialization...");
        // make sure planner initialized
        auto traj_ptr = traj_port->get_new<AccTrajectory>();
        t0 = std::chrono::system_clock::now();
        auto state_ptr = state_port->get_new<Measurement>();
    }

    /// loop
    logger.info("Controller Block {} started", get_block_name());
    while (true)
    {
        // retrieving values
        auto state = state_port->get_new<Measurement>();

        VectorXd pdes, vdes, ades;
        auto traj = get_desired(pdes, vdes, ades, state->timestamp);

        // time between controller call
        auto now = std::chrono::system_clock::now();
        double dt =
            std::chrono::duration_cast<std::chrono::microseconds>(now - t_prev)
                .count() /
            1e6;
        t_prev = now;

        // update errors
        VectorXd e = pdes - state->pos;
        VectorXd ed = vdes - state->vel;

        for (uint32_t i = 0; i < rob_model.num_dof; i++)
        {
            while (e[i] < -M_PI)
                e[i] += 2 * M_PI;
            while (e[i] > M_PI)
                e[i] -= 2 * M_PI;
        }

        // call for dynamics calculation
        InputPack inputs{pdes,       vdes, traj->acc, state->pos,
                         state->vel, e,    ed,        dt};
        ctrl.control = dynamics(inputs);
        ctrl.is_gripper_open = traj->is_gripper_open;

        // update output
        ctrl.frame_id = state->frame_id + 1;
        default_output->set<ControlMsg>(ctrl);
    }
}

ControllerBlock::ControllerBlock(const std::string &block_name,
                                 const std::string &model_path)
    : ControllerBlock(block_name, Model::Parser::parse(model_path))
{
}

ControllerBlock::ControllerBlock(const std::string &block_name,
                                 const Model::SpatialModel<double> &rob_model)
    : ::Roahm::System::BaseBlock(block_name,
                                 System::OutputPort<ControlMsg>::make_shared()),
      traj_port(register_input<AccTrajectory>("trajectory")),
      state_port(register_input<Measurement>("state")), rob_model(rob_model),
      t0(std::chrono::system_clock::now())
{
}

std::shared_ptr<const msgs::AccTrajectory> ControllerBlock::get_desired(
    VectorXd &pdes, VectorXd &vdes, VectorXd &ades,
    decltype(std::chrono::system_clock::now()) t)
{
    using namespace std::chrono_literals;
    // check new trajectory
    auto latest_traj = traj_port->get_nowait<AccTrajectory>();
    // add to queue if updated
    if (!has_trajectory || latest_traj.get() != trajectories.back().get())
    {
        // initialization
        if (!has_trajectory || latest_traj->reset)
        {
            // clear all trajectory
            trajectories.clear();
            has_trajectory = true;
            t0 = t;
            // NOTE: deal with condition when start_time!=0
            std::size_t duration =
                std::chrono::time_point_cast<std::chrono::milliseconds>(t0)
                    .time_since_epoch()
                    .count();

            std::chrono::milliseconds dur(
                duration - std::size_t(latest_traj->start_time * 1000));

            t0 = std::chrono::time_point<std::chrono::system_clock>(dur);
        }
        trajectories.push_back(latest_traj);
    }

    // check validity of the first trajectory
    {
        const double t_rel =
            std::chrono::duration_cast<std::chrono::microseconds>(t - t0)
                .count() /
            1e6; // time passed

        while (!trajectories.empty())
        {
            auto first_traj = trajectories.front();
            // not yet reached? NOTE: this should only happen if planner has
            // error
            if (t_rel < first_traj->start_time)
            {
                throw std::runtime_error(
                    "Planner Error: Oldest trajectory still in the future!");
            }
            // expired?
            else if (t_rel > (first_traj->duration + first_traj->start_time))
            {
                trajectories.pop_front();
            }
            // passed the check!
            else
            {
                break;
            }
        }
    }

    // fetch valid trajectory
    if (trajectories.empty())
    {
        throw std::runtime_error("Planner Error: no new trajectories!");
    }

    // get desired trajectory
    auto traj = trajectories.front();
    const double t_rel =
        std::chrono::duration_cast<std::chrono::microseconds>(t - t0).count() /
            1e6 -
        traj->start_time; // time passes since start time

    // handle based on traj type difference
    if (traj->is_bernstein)
    {
        pdes = Helper::bernstein_q_des(traj->pos, traj->vel, traj->acc, traj->k,
                                       t_rel, traj->trajectory_duration);
        vdes =
            Helper::bernstein_qd_des(traj->pos, traj->vel, traj->acc, traj->k,
                                     t_rel, traj->trajectory_duration);
        ades =
            Helper::bernstein_qdd_des(traj->pos, traj->vel, traj->acc, traj->k,
                                      t_rel, traj->trajectory_duration);
    }
    else
    {
        ades = traj->acc;
        vdes = ades * t_rel + traj->vel;
        pdes = (vdes + traj->vel) * t_rel / 2 + traj->pos;
    }

    return traj;
}


void ControllerBlock::set_trajectories(
    std::deque<std::shared_ptr<const AccTrajectory>> trajs)
{
    trajectories = trajs;
}

std::shared_ptr<const msgs::AccTrajectory> ControllerBlock::get_desired_csv(
    VectorXd &pdes, VectorXd &vdes, VectorXd &ades,
    decltype(std::chrono::system_clock::now()) t)
{
    using namespace std::chrono_literals;

    // check validity of the first trajectory
    {
        const double t_rel =
            std::chrono::duration_cast<std::chrono::microseconds>(t - t0)
                .count() /
            1e6; // time passed

        while (!trajectories.empty())
        {
            auto first_traj = trajectories.front();
            // not yet reached? NOTE: this should only happen if planner has
            // error
            if (t_rel < first_traj->start_time)
            {
                throw std::runtime_error(
                    "Planner Error: Oldest trajectory still in the future!");
            }

            // expired?
            else if (t_rel > (first_traj->duration + first_traj->start_time))
            {
                trajectories.pop_front();
            }

            // passed the check!
            else
            {
                break;
            }
        }
    }

    // fetch valid trajectory
    if (trajectories.empty())
    {
        throw std::runtime_error("Planner Error: no new trajectories!");
    }

    // get desired trajectory
    auto traj = trajectories.front();
    const double t_rel =
        std::chrono::duration_cast<std::chrono::microseconds>(t - t0).count() /
            1e6 -
        traj->start_time; // time passes since start time

    // handle based on traj type difference
    if (traj->is_bernstein)
    {
        pdes = Helper::bernstein_q_des(traj->pos, traj->vel, traj->acc, traj->k,
                                       t_rel, traj->trajectory_duration);
        vdes =
            Helper::bernstein_qd_des(traj->pos, traj->vel, traj->acc, traj->k,
                                     t_rel, traj->trajectory_duration);
        ades =
            Helper::bernstein_qdd_des(traj->pos, traj->vel, traj->acc, traj->k,
                                      t_rel, traj->trajectory_duration);
    }
    else
    {
        ades = traj->acc;
        vdes = ades * t_rel + traj->vel;
        pdes = (vdes + traj->vel) * t_rel / 2 + traj->pos;
    }

    return traj;

}
} // namespace Dynamics
} // namespace Roahm
