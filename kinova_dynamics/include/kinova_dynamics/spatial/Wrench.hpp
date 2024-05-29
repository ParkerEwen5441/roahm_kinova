/**
 * @brief declaration for force spatial vector
 * @file
 **/
#pragma once
#include "SpatialBase.hpp"

namespace Roahm
{
namespace Spatial
{
/**
 * @brief Wrench spatial vector class
 **/
template <typename T> class Wrench : public SpatialBase<T, Wrench<T>>
{
  protected:
    using Base = SpatialBase<T, Wrench<T>>;

  public:
    // parent constructor
    using Base::Base;
};
} // namespace Spatial
} // namespace Roahm
