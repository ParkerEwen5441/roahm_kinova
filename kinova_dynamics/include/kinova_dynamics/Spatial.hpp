#pragma once
// includes
#include "spatial/Interval.hpp"
#include "spatial/SpatialAxis.hpp"
#include "spatial/SpatialInertia.hpp"
#include "spatial/SpatialTransform.hpp"
#include "spatial/Twist.hpp"
#include "spatial/Wrench.hpp"
#include <cmath>

namespace Roahm
{
namespace Spatial
{
/** Typedef for convenience */
// double
using Twistd = Twist<double>;
using Wrenchd = Wrench<double>;
using SAxisd = SpatialAxis<double>;
using SInertiad = SpatialInertia<double>;
using STransformd = SpatialTransform<double>;

// interval
using Intd = Interval<double>;
using IntTwistd = Twist<Intd>;
using IntWrenchd = Wrench<Intd>;
using IntSInertiad = SpatialInertia<Intd>;
using IntSTransformd = SpatialTransform<Intd>;
using IntSAxisd = SpatialAxis<Intd>;

/**
 * @brief implement math functions
 **/
Interval<double> sin(Interval<double> val);
Interval<double> cos(Interval<double> val);

/**
 * @brief output a rotation matrix based on directional vector
 * @param v_dir: direction vector
 **/
template <typename T>
Eigen::Matrix<T, 3, 3> dir2rot(T angle, const Eigen::Matrix<T, 3, 1> &axis)
{
    using Matrix3 = Eigen::Matrix<T, 3, 3>;
    using Vector3 = Eigen::Matrix<T, 3, 1>;
    using std::cos, std::sin;
    Matrix3 res;
    Vector3 sin_axis = sin(angle) * axis;
    T c = cos(angle);
    Vector3 cos1_axis = (T(1) - c) * axis;

    T tmp;
    tmp = cos1_axis.x() * axis.y();
    res.coeffRef(0, 1) = tmp - sin_axis.z();
    res.coeffRef(1, 0) = tmp + sin_axis.z();

    tmp = cos1_axis.x() * axis.z();
    res.coeffRef(0, 2) = tmp + sin_axis.y();
    res.coeffRef(2, 0) = tmp - sin_axis.y();

    tmp = cos1_axis.y() * axis.z();
    res.coeffRef(1, 2) = tmp - sin_axis.x();
    res.coeffRef(2, 1) = tmp + sin_axis.x();

    res.diagonal() = (cos1_axis.cwiseProduct(axis)).array() + c;
    return res;
}
} // namespace Spatial
} // namespace Roahm
