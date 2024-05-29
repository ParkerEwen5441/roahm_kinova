#pragma once
#include "kinova_dynamics/Model.hpp"
#include "kinova_dynamics/Parser.hpp"
#include "kinova_dynamics/Spatial.hpp"
#include "kinova_dynamics/spatial/SpatialAxis.hpp"
#include "kinova_dynamics/spatial/SpatialTransform.hpp"
#include <eigen3/Eigen/Core>

namespace Roahm
{
namespace Dynamics
{
using ::Roahm::Model::SpatialLink;
using namespace Spatial;

/**
 * @brief determine friction direction and apply friction
 * @param friction: friction value to apply
 * @param v: joint velocity
 * @param a: joint acceleration
 **/
template <typename T>
T apply_friction(T friction, double v, double a, double thr = 0.05)
{
    // static case
    if (std::abs(v) < thr)
    {
        // no motion, no friction
        if (std::abs(a) < thr)
            return 0;
        return friction * (a > 0 ? 1. : -1.);
    }

    // non-static, just oppsite to moving direction
    else
    {
        return friction * (v > 0 ? 1. : -1.);
    }
}

/**
 * @brief calculate rnea
 * @param robot: robot model
 * @param q: joint position
 * @param qd: joint velocity
 * @param qdd: joint acceleration
 * @return calculated torque
 **/
template <typename T>
VecX<T> rnea_T(const Model::SpatialModel<T> &robot, const VectorXd &q,
               const VectorXd &qd, const VectorXd &qd_aux, const VectorXd &qdd,
               bool use_gravity = true)
{
    // namespace setting
    using JointType = Model::JointType;
    using namespace ::Roahm::Model;
    using namespace ::Roahm::Spatial;
    using std::vector;
    // alias
    auto num_dof = robot.num_dof;
    const std::vector<SpatialLink<T>> &links = robot.links;
    // check
    assert(q.size() == num_dof);
    assert(qd.size() == num_dof);
    assert(qdd.size() == num_dof);
    vector<Twist<T>> v(links.size());
    vector<Twist<T>> v_aux(links.size());
    vector<Twist<T>> a(links.size());
    vector<Wrench<T>> f(links.size());
    vector<SpatialTransform<T>> Xs_li_i(links.size());
    
    // spatial axis w.r.t. base
    vector<SpatialAxis<T>> Sb;
    Sb.reserve(links.size());
    VecX<T> tau(num_dof);

    /** Forward pass */
    uint32_t di = 0; // index for dof
    for (uint32_t i = 0; i < links.size(); i++)
    {
        // get parent index
        const SpatialLink<T> &link = links[i];
        // NOTE: cope with previous link
        Twist<T> pv{}, pa{}, pv_aux{};
        if (i != 0)
        {
            pv = v[links[i].parent_id];
            pv_aux = v_aux[links[i].parent_id];
            pa = a[links[i].parent_id];
        }
        else if (use_gravity)
        {
            pa.linear = -robot.gravity;
        }

        // cope with fixed type
        double jq{0}, jqd{0}, jqdd{0}, jqd_aux{0};
        if (link.type != JointType::FIXED)
        {
            jq = q[di];
            jqd = qd[di];
            jqd_aux = qd_aux[di];
            jqdd = qdd[di];
            di++;
        }

        // joint axis
        Twist<T> s = link.Xbw.inv_apply(link.S.as_twist());
        Sb.push_back(SpatialAxis<T>(s.angular, s.linear));

        // joint transform
        SpatialTransform<T> &Xli_i = Xs_li_i[i];
        Xli_i = link.get_transform(jq, Sb[i]);

        /** update kinematics */
        // NOTE: for baselink, v[li] = v[i] = 0
        // v = previous + joint movement
        Twist<T> jvel = Sb[i] * jqd;
        Twist<T> jvel_aux = Sb[i] * jqd_aux;

        // v = previous + joint
        v[i] = Xli_i.apply(pv) + jvel;
        v_aux[i] = Xli_i.apply(pv_aux) + jvel_aux;

        // a = previous + joint acceleration + coriolis term
        a[i] = Xli_i.apply(pa) + Sb[i] * jqdd + v[i].cross(jvel_aux);

        // calculate v x Iv Jon's way
        Wrench<T> vIv{v_aux[i].angular.cross(link.SI.I() * v[i].angular) +
                          link.SI.I() * v_aux[i].angular.cross(v[i].angular),
                      link.SI.mass() * v_aux[i].angular.cross(v[i].linear)};

        // // force update
        f[i] = link.SI * a[i] + vIv;
    }

    /** Backward pass */
    for (uint32_t i = links.size() - 1, di = num_dof - 1; i < links.size(); i--)
    {
        auto &link = links[i];
        if (link.type != JointType::FIXED)
        {
            tau[di] = Sb[i].dot(f[i]);
            // damping
            tau[di] += link.damping * qd[di];
            // friction
            tau[di] += apply_friction(link.friction, qd[di], qdd[di]);
            // amature
            tau[di] += link.armature * qdd[di];
            di--;
        }
        if (i != 0)
        {
            f[link.parent_id] += Xs_li_i[i].T_apply(f[i]);
        }
    }
    return tau;
}

/**
 * @brief calculate mass matrix using RNEA
 * @param robot: robot model to be used
 * @param q: robot configuration
 * @return mass matrix of shape (num_dof, num_dof)
 **/
template <typename T>
MatX<T> mass_matrix(const Model::SpatialModel<T> &robot, const VectorXd &q)
{
    // check dimension
    assert(q.size() == robot.num_dof);
    MatX<T> mass_mat(robot.num_dof, robot.num_dof);
    const VectorXd zero_state = VectorXd::Zero(robot.num_dof);
    VectorXd qdd(robot.num_dof);
    for (std::size_t i = 0; i < robot.num_dof; i++)
    {
        qdd.setZero();
        qdd(i) = 1;
        mass_mat.block(0, i, robot.num_dof, 1) =
            rnea_T(robot, q, zero_state, zero_state, qdd, false);
    }
    return mass_mat;
}

/**
 * @brief wrapper for normal rnea
 **/
inline auto rnea(const Model::SpatialModel<double> &robot, const VectorXd &q,
                 const VectorXd &qd, const VectorXd &qd_aux,
                 const VectorXd &qdd, bool use_gravity = true)
{
    return rnea_T<double>(robot, q, qd, qd_aux, qdd, use_gravity);
}

inline auto int_rnea(const Model::SpatialModel<Intd> &robot, const VectorXd &q,
                     const VectorXd &qd, const VectorXd &qd_aux,
                     const VectorXd &qdd, bool use_gravity = true)
{
    return rnea_T<Intd>(robot, q, qd, qd_aux, qdd, use_gravity);
}
} // namespace Dynamics
} // namespace Roahm
