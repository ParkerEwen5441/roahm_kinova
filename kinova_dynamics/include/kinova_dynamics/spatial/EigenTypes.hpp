/**
 * @brief custom eigen typedefs
 * @file
 **/
#pragma once
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/src/Core/util/Constants.h>

namespace Roahm
{
namespace Spatial
{
template <typename T> using Vec2 = Eigen::Matrix<T, 2, 1>;
template <typename T> using Vec3 = Eigen::Matrix<T, 3, 1>;
template <typename T> using Vec6 = Eigen::Matrix<T, 6, 1>;
template <typename T> using VecX = Eigen::Matrix<T, Eigen::Dynamic, 1>;
template <typename T> using Mat3 = Eigen::Matrix<T, 3, 3>;
template <typename T>
using MatX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
} // namespace Spatial
} // namespace Roahm
