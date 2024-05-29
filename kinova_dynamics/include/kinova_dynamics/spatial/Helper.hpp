/**
 * @file
 * @store helper functions for spatial lib
 **/
#pragma once
#include "EigenTypes.hpp"
#include <iostream>

namespace Roahm
{
namespace Spatial
{
/**
 * @brief get skew symmetric matrix representation of vector
 * @return 3x3 matrix
 **/
template <typename T> Mat3<T> skew_sym_mat(const Vec3<T> &v)
{
    Mat3<T> v_hat;
    v_hat << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
    return v_hat;
}

} // namespace Spatial
} // namespace Roahm
