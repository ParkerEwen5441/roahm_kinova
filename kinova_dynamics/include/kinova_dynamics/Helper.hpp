#pragma once
#include <eigen3/Eigen/Eigen>
namespace Roahm
{
namespace Dynamics
{
namespace Helper
{
using Eigen::VectorXd;

/**
 * @brief calculate q for bernstein traj
 * @param q0: initial q
 * @param qd0: initial qd
 * @param qdd0: initial qdd
 * @param t: time
 * @param traj_duration: duration of the trajectory, default 1
 **/
VectorXd bernstein_q_des(const VectorXd &q0, VectorXd qd0, VectorXd qdd0,
                         const VectorXd &k, double t,
                         const double traj_duration = 1.0);

/**
 * @brief calculate qd for bernstein traj
 * @param q0: initial q
 * @param qd0: initial qd
 * @param qdd0: initial qdd
 * @param t: time
 * @param traj_duration: duration of the trajectory, default 1
 **/
VectorXd bernstein_qd_des(const VectorXd &q0, VectorXd qd0, VectorXd qdd0,
                          const VectorXd &k, double t,
                          const double traj_duration = 1.0);

/**
 * @brief calculate qdd for bernstein traj
 * @param q0: initial q
 * @param qd0: initial qd
 * @param qdd0: initial qdd
 * @param t: time
 * @param traj_duration: duration of the trajectory, default 1
 **/
VectorXd bernstein_qdd_des(const VectorXd &q0, VectorXd qd0, VectorXd qdd0,
                           const VectorXd &k, double t,
                           const double traj_duration = 1.0);
} // namespace Helper
} // namespace Dynamics
} // namespace Roahm
