/**
 * @file
 * @brief Forward Dynamics through Composite Rigid Body Algorithm
 **/
#pragma once
#include "kinova_dynamics/RNEA.hpp"
#include <eigen3/Eigen/src/Core/util/Constants.h>

namespace Roahm
{
namespace Dynamics
{
/**
 * @brief calculate forward dynamics using CRBA
 * NOTE: Work In Progress
 **/
template <typename T>
VectorXd crba_T(const Model::SpatialModel<T> &robot, const VectorXd &q,
                const VectorXd &qd, const VecX<T> &torque)
{
    // alias for easy reference
    auto num_dof = robot.num_dof;
    /** C */
    VectorXd qdd = VectorXd(num_dof);
    qdd.setZero();
    auto C = rnea_T<T>(robot, q, qd, qd, qdd);

    /** H */
    std::size_t num_links = robot.links.size();
    auto &links = robot.links;
    std::vector<STransformd> Xs_li_i(links.size());
    uint32_t di = 0; // index for dof
    for (uint32_t i = 0; i < links.size(); i++)
    {
        // cope with fixed type
        double jq{0};
        if (links[i].type != Model::JointType::FIXED)
        {
            jq = q[di];
            di++;
        }

        // joint transform
        Xs_li_i[i] = links[i].get_transform(jq);
    }
    // composite inertia
    std::vector<SpatialInertia<T>> IC;
    for (const auto &link : robot.links)
    {
        IC.push_back(link.SI);
    }

    for (int i = num_links; i >= 0; i--)
    {
        if (i != 0)
        {
            IC[robot.links[i].parent_id] += IC[i].transform(Xs_li_i[i]);
        }
    }

    // H
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> H;
    H.resize(num_links, num_links);
    // for (std::size_t i = 0; i < num_links; i++)
    // {
    //     IC[i] * links[i].S
    // }
}
} // namespace Dynamics
} // namespace Roahm
