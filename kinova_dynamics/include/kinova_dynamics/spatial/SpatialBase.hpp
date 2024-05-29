#pragma once
#include "EigenTypes.hpp"
#include <iostream>

namespace Roahm
{
namespace Spatial
{
using namespace Eigen;

/**
 * @brief base template class for spatial vector, can't be initialized directly
 * @details using CRTP to overload functions
 * @param T: basic arithmetic type for eigen vectors to hold
 **/
template <typename T, typename Derived> class SpatialBase
{
    // methods
  public:
    /** Constructors */
    /**
     * @brief default constructor
     * @details perform setZero()
     **/
    SpatialBase()
    {
        angular.setZero();
        linear.setZero();
    }

    /**
     * @brief construct from angular and linear part
     **/
    SpatialBase(const Vec3<T> &w, const Vec3<T> &v) : angular{w}, linear{v}
    {
    }

    SpatialBase(Vec3<T> &&w, Vec3<T> &&v)
        : angular{std::move(w)}, linear{std::move(v)}
    {
    }

    /**
     * @brief constructor from rvalue
     **/
    SpatialBase(Derived &&sv)
        : angular{std::move(sv.angular)}, linear{std::move(sv.linear)}
    {
    }

    /** OPERATORS */
    /**
     * @brief overload for -V
     * @return new instance of {-angular, -linear}
     **/
    Derived operator-()
    {
        return Derived{-angular, -linear};
    }

    /**
     * @brief overload for = using copy
     * @param sv: operator to copy
     **/
    Derived &operator=(const Derived &sv)
    {
        if (this != &sv)
        {
            angular = sv.angular;
            linear = sv.linear;
        }
        return *static_cast<Derived *>(this);
    }

    /**
     * @brief overload for = using move
     * @param sv: vector to move
     **/
    Derived &operator=(Derived &&sv)
    {
        if (this != &sv)
        {
            angular = std::move(sv.angular);
            linear = std::move(sv.linear);
        }
        return *static_cast<Derived *>(this);
    }

    /**
     * @brief overload +
     * @details add angular and linear part separately
     * @param sv: another spatial vector
     * @return new instance of sum of two %SpatialBase
     **/
    Derived operator+(const Derived &sv) const
    {
        return Derived{angular + sv.angular, linear + sv.linear};
    }

    /**
     * @brief overload +=
     * @details perform inplace addition
     * @param sv: another spatial vector
     * @return reference to the instance
     **/
    Derived &operator+=(const Derived &sv)
    {
        angular += sv.angular;
        linear += sv.linear;
        return *static_cast<Derived *>(this);
    }

    /**
     * @brief overload -
     * @details substract angular and linear part separately
     * @param sv: another spatial vector
     * @return new instance of substraction of two %SpatialBase
     **/
    Derived operator-(const Derived &sv) const
    {
        return Derived{angular - sv.angular, linear - sv.linear};
    }

    /**
     * @brief overload -=
     * @details perform inplace substraction
     * @param sv: another spatial vector
     * @return reference to the instance
     **/
    Derived &operator-=(const Derived &sv)
    {
        angular -= sv.angular;
        linear -= sv.linear;
        return *static_cast<Derived *>(this);
    }

    /**
     * @brief overload for * with scalar
     * @param t: scalar to be multiplied
     * @return new Derived class
     **/
    Derived operator*(const T &t) const
    {
        return Derived{angular * t, linear * t};
    }

    friend Derived operator*(const T &t, const Derived &sv)
    {
        return sv * t;
    }

    /**
     * @brief overload for * with scalar
     * @param t: scalar to be multiplied
     * @return *this
     **/
    Derived &operator*=(const T &t)
    {
        angular *= t;
        linear *= t;
        return *static_cast<Derived *>(this);
    }

    /**
     * @brief overload for / with scalar
     * @param t: scalar to divide
     * @return new Derived class
     **/
    Derived operator/(const T &t) const
    {
        return Derived{angular / t, linear / t};
    }

    /**
     * @brief overload for /= with scalar
     * @param t: scalar to divide
     * @return *this
     **/
    Derived &operator/=(const T &t)
    {
        angular /= t;
        linear /= t;
        return *static_cast<Derived *>(this);
    }

    /**
     * @brief overload for ==
     * @details compare separately
     * @return true if both part match
     **/
    bool operator==(const Derived &sv) const
    {
        return (angular == sv.angular) && (linear == sv.linear);
    }

    /**
     * @brief overload for !=
     * @details compare separately
     * @return true if both part match
     **/
    bool operator!=(const Derived &sv) const
    {
        return (angular != sv.angular) || (linear != sv.linear);
    }

    /** OTHER OPERATIONS */
    /**
     * @brief dot product with another %SpatialBase
     * @param sv: another spatial vector
     * @return result in basic type T
     **/
    T dot(const Derived &sv) const
    {
        return angular.dot(sv.angular) + linear.dot(sv.linear);
    }

    // variables
  public:
    Vec3<T> angular;
    Vec3<T> linear;
};

/**
 * @brief ostream for %SpatialBase
 **/
template <typename T, typename Derived>
std::ostream &operator<<(std::ostream &os, const SpatialBase<T, Derived> &sv)
{
    os << "angular: " << sv.angular.transpose() << '\n';
    os << "linear: " << sv.linear.transpose();
    return os;
}
} // namespace Spatial
} // namespace Roahm
