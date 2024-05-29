#include "kinova_dynamics/Helper.hpp"

namespace Roahm
{
namespace Dynamics
{
namespace Helper
{
VectorXd bernstein_q_des(const VectorXd &q0, VectorXd qd0, VectorXd qdd0,
                         const VectorXd &k, double t,
                         const double traj_duration)
{
    qd0 = qd0 * traj_duration;
    qdd0 = qdd0 * traj_duration * traj_duration;
    t = t / traj_duration;

    double B0 = -pow(t - 1, 5);
    double B1 = 5 * t * pow(t - 1, 4);
    double B2 = -10 * pow(t, 2) * pow(t - 1, 3);
    double B3 = 10 * pow(t, 3) * pow(t - 1, 2);
    double B4 = -5 * pow(t, 4) * (t - 1);
    double B5 = pow(t, 5);
    VectorXd beta0 = q0;
    VectorXd beta1 = q0 + qd0 / 5;
    VectorXd beta2 = q0 + (2 * qd0) / 5 + qdd0 / 20;
    VectorXd beta3 = q0 + k;
    VectorXd beta4 = q0 + k;
    VectorXd beta5 = q0 + k;
    return B0 * beta0 + B1 * beta1 + B2 * beta2 + B3 * beta3 + B4 * beta4 +
           B5 * beta5;
}

VectorXd bernstein_qd_des(const VectorXd &q0, VectorXd qd0, VectorXd qdd0,
                          const VectorXd &k, double t,
                          const double traj_duration)
{
    qd0 = qd0 * traj_duration;
    qdd0 = qdd0 * traj_duration * traj_duration;
    t = t / traj_duration;

    double dB0 = pow(t - 1.0, 4.0) * -5.0;
    double dB1 = t * pow(t - 1.0, 3.0) * 2.0E+1 + pow(t - 1.0, 4.0) * 5.0;
    double dB2 =
        t * pow(t - 1.0, 3.0) * -2.0E+1 - (t * t) * pow(t - 1.0, 2.0) * 3.0E+1;
    double dB3 = pow(t, 3.0) * (t * 2.0 - 2.0) * 1.0E+1 +
                 (t * t) * pow(t - 1.0, 2.0) * 3.0E+1;
    double dB4 = pow(t, 3.0) * (t - 1.0) * -2.0E+1 - pow(t, 4.0) * 5.0;
    double dB5 = pow(t, 4.0) * 5.0;
    VectorXd beta0 = q0;
    VectorXd beta1 = q0 + qd0 / 5;
    VectorXd beta2 = q0 + (2 * qd0) / 5 + qdd0 / 20;
    VectorXd beta3 = q0 + k;
    VectorXd beta4 = q0 + k;
    VectorXd beta5 = q0 + k;
    return (dB0 * beta0 + dB1 * beta1 + dB2 * beta2 + dB3 * beta3 +
            dB4 * beta4 + dB5 * beta5) /
           traj_duration;
}

VectorXd bernstein_qdd_des(const VectorXd &q0, VectorXd qd0, VectorXd qdd0,
                           const VectorXd &k, double t,
                           const double traj_duration)
{
    qd0 = qd0 * traj_duration;
    qdd0 = qdd0 * traj_duration * traj_duration;
    t = t / traj_duration;

    double t2 = t * 2.0;
    double t3 = t * t;
    double t4 = t * t * t;
    double t5 = t - 1.0;
    double t6 = t2 - 2.0;
    double t7 = t4 * 2.0E+1;
    double t8 = t5 * t5;
    double t9 = t5 * t5 * t5;
    double t10 = t9 * 2.0E+1;
    double t11 = t * t8 * 6.0E+1;
    double t12 = -t10;
    double ddB0 = t12;
    double ddB1 = t9 * 4.0E+1 + t11;
    double ddB2 = t12 - t * t8 * 1.2E+2 - t3 * t6 * 3.0E+1;
    double ddB3 = t7 + t11 + t3 * t6 * 6.0E+1;
    double ddB4 = t4 * -4.0E+1 - t3 * t5 * 6.0E+1;
    double ddB5 = t7;
    VectorXd beta0 = q0;
    VectorXd beta1 = q0 + qd0 / 5;
    VectorXd beta2 = q0 + (2 * qd0) / 5 + qdd0 / 20;
    VectorXd beta3 = q0 + k;
    VectorXd beta4 = q0 + k;
    VectorXd beta5 = q0 + k;
    return (ddB0 * beta0 + ddB1 * beta1 + ddB2 * beta2 + ddB3 * beta3 +
            ddB4 * beta4 + ddB5 * beta5) /
           (traj_duration * traj_duration);
}
} // namespace Helper
} // namespace Dynamics
} // namespace Roahm
