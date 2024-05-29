/**
 * @brief motion spatial vector
 * @file
 **/
#pragma once
#include "SpatialBase.hpp"
#include "Wrench.hpp"

namespace Roahm
{
namespace Spatial
{
/**
 * @brief motion spatial vector
 * @param T: basic value stored in matrix
 **/
template <typename T> class Twist : public SpatialBase<T, Twist<T>>
{
  private:
    using Base = SpatialBase<T, Twist<T>>;

  public:
    // inherit constructors
    using Base::Base;

    /**
     * @brief cross product with another Twist
     **/
    Twist<T> cross(const Twist<T> &t) const
    {
        return Twist{this->angular.cross(t.angular),
                     this->angular.cross(t.linear) +
                         this->linear.cross(t.angular)};
    }

    /**
     * @brief cross product with %Wrench
     * @param w: %Wrench to cross
     * @return new %Wrench
     **/
    Wrench<T> cross_w(const Wrench<T> &w) const
    {
        return Wrench<T>{this->angular.cross(w.angular) +
                             this->linear.cross(w.linear),
                         this->angular.cross(w.linear)};
    }
};
} // namespace Spatial
} // namespace Roahm
