#pragma once
#include "SpatialBase.hpp"
#include "Twist.hpp"
#include "Wrench.hpp"
#include <stdexcept>

namespace Roahm
{
namespace Spatial
{
template <typename T> class SpatialAxis : private SpatialBase<T, SpatialAxis<T>>
{
  private:
    using Base = SpatialBase<T, SpatialAxis<T>>;

  public:
    SpatialAxis() : Base()
    {
    }

    /**
     * @brief Ctor indicating index of axis */
    SpatialAxis(uint32_t index) : SpatialAxis()
    {
        // angular
        if (index < 3)
        {
            this->angular[index] = 1;
        }
        // linear
        else if (index < 6)
        {
            this->linear[index - 3] = 1;
        }
        else
        {
            throw std::out_of_range(
                "index out of spatial vector with dimension 6");
        }
    }

    /**
     * @brief constructor using only one
     **/
    SpatialAxis(const Vec3<T> &v, bool is_linear) : SpatialAxis()
    {
        if (is_linear)
            this->linear = v.normalized();
        else
            this->angular = v.normalized();
    }

    SpatialAxis(const Vec3<T> &angular, const Vec3<T> &linear)
        : Base(angular, linear)
    {
    }

    /**
     * @brief copy Ctor
     **/
    SpatialAxis(const SpatialAxis<T> &s) : Base(s.angular, s.linear)
    {
    }

    /**
     * @brief treat as twist to enable more operations
     **/
    Twist<T> as_twist() const
    {
        return Twist<T>{this->angular, this->linear};
    }

    /**
     * @brief operator=
     **/
    SpatialAxis<T> &operator=(const SpatialAxis<T> s)
    {
        this->angular = s.angular;
        this->linear = s.linear;
        return *this;
    }

    /**
     * @brief constant angulalr part
     **/
    const Vec3<T> &ang() const
    {
        return this->angular;
    }

    /**
     * @brief constant linear part
     **/
    const Vec3<T> &lin() const
    {
        return this->linear;
    }

    /**
     * @brief get twist by multiply a scalar
     **/
    Twist<T> operator*(T val) const
    {
        return Twist<T>(this->angular * val, this->linear * val);
    }

    /**
     * @brief dot product with wrench
     **/
    T dot(const Wrench<T> &w) const
    {
        return this->linear.dot(w.linear) + this->angular.dot(w.angular);
    }

    /**
     * @brief overload for output
     **/
    friend std::ostream &operator<<(std::ostream &os, const SpatialAxis<T> &sa)
    {
        os << "angular: " << sa.angular.transpose() << '\n';
        os << "linear: " << sa.linear.transpose();
        return os;
    }

    /**
     * @brief Ctor from twist
     **/
    // SpatialAxis()

    /**
     * @brief get matrix that maps vector onto the axis
     **/
    // Eigen::Matrix<> mat_form() const
    // {
    //     Mat3<T> ret;
    //     ret <<
    // }
};
} // namespace Spatial
} // namespace Roahm
