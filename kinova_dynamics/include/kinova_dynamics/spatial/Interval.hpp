/**
 * @brief Implementation of Interval Scalar type
 * @file
 * @ref https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node45.html
 **/
#pragma once
#include <algorithm>
#include <chrono>
#include <cmath>
#include <eigen3/Eigen/Geometry>
#include <limits>
#include <ostream>
#include <utility>

namespace Roahm
{
namespace Spatial
{

/**
 * @brief interval type
 * @param T: basic scalar type
 **/
template <typename T> class Interval
{
  public:
    // lower and upper bound
    T lower, upper;

  public:
    /** Constructors */
    /**
     * @brief default constructor
     * @details init lower and upper by their default value
     **/
    Interval() : lower{}, upper{}
    {
    }

    /**
     * @brief constructor from single scalar
     **/
    Interval(T val) : lower{val}, upper{val}
    {
    }

    /**
     * @brief constructor
     * @param lower: lower bound
     * @param upper: upper bound
     **/
    Interval(T lower, T upper) : lower{lower}, upper{upper}
    {
    }

    /**
     * @brief copy constructor
     * @param in: another interval
     **/
    Interval(const Interval<T> &in) : lower{in.lower}, upper{in.upper}
    {
    }

    /**
     * @brief copy constructor from rvalue ref
     * @param in: another interval
     **/
    Interval(Interval<T> &&in)
        : lower{std::move(in.lower)}, upper{std::move(in.upper)}
    {
    }

    /** Operators */
    // 1
    /**
     * @brief -int
     * @return -int
     **/
    Interval<T> operator-() const
    {
        return Interval<T>{-upper, -lower};
    }

    // 2
    /**
     * @brief operator =
     * @param in: another interval
     * @return *this
     **/
    Interval<T> &operator=(const Interval<T> &in)
    {
        lower = in.lower;
        upper = in.upper;
        return *this;
    }

    Interval<T> &operator=(Interval<T> &&in)
    {
        lower = std::move(in.lower);
        upper = std::move(in.upper);
        return *this;
    }

    /**
     * @brief +
     * @param in: another interval
     * @return int1 + int2
     **/
    Interval<T> operator+(const Interval<T> &in) const
    {
        return Interval<T>{lower + in.lower, upper + in.upper};
    }

    /**
     * @brief +=
     * @param in: another interval
     * @return *this
     **/
    Interval<T> &operator+=(const Interval<T> &in)
    {
        lower += in.lower;
        upper += in.upper;
        return *this;
    }

    /**
     * @brief -
     * @param in: another interval
     * @return int1 - int2
     **/
    Interval<T> operator-(const Interval<T> &in) const
    {
        return Interval<T>{lower - in.upper, upper - in.lower};
    }

    /**
     * @brief -=
     * @param in: another interval
     * @return *this
     **/
    Interval<T> &operator-=(const Interval<T> &in)
    {
        // NOTE: deal with &in == this
        this->operator=((*this) - in);
        return *this;
    }

    /**
     * @brief *
     * @param in: another interval
     * @return int1 * int2
     **/
    Interval<T> operator*(const Interval<T> &in) const
    {
        return Interval<T>{
            std::min(std::min(lower * in.lower, lower * in.upper),
                     std::min(upper * in.lower, upper * in.upper)),
            std::max(std::max(lower * in.lower, lower * in.upper),
                     std::max(upper * in.lower, upper * in.upper))};
    }

    /**
     * @brief *=
     * @param in: another interval
     * @return *this
     **/
    Interval<T> &operator*=(const Interval<T> &in)
    {
        // assign rvalue
        this->operator=((*this) * in);
        return *this;
    }

    /**
     * @brief /
     * @param in: another interval
     * @return int1 / int2
     **/
    Interval<T> operator/(const Interval<T> &in) const
    {
        return Interval<T>{
            std::min(std::min(lower / in.lower, lower / in.upper),
                     std::min(upper / in.lower, upper / in.upper)),
            std::max(std::max(lower / in.lower, lower / in.upper),
                     std::max(upper / in.lower, upper / in.upper))};
    }

    Interval<T> &operator/=(const Interval<T> &in)
    {
        // assign rvalue
        this->operator=((*this) / in);
        return *this;
    }

    /** Operation with scalar */
    /**
     * @brief * with scalar
     **/
    Interval<T> operator*(const T &t)
    {
        return Interval<T>{t * lower, t * upper};
    }

    friend inline Interval<T> operator*(const T &t, const Interval<T> &in)
    {
        return in * t;
    }

    /**
     * @brief *= with scalar
     **/
    Interval<T> &operator*=(const T &t)
    {
        lower *= t;
        upper *= t;
        return *this;
    }

    /**
     * @brief / with scalar
     **/
    Interval<T> operator/(const T &t)
    {
        return Interval<T>{lower / t, upper / t};
    }

    /**
     * @brief /= with scalar
     **/
    Interval<T> &operator/=(const T &t)
    {
        lower /= t;
        upper /= t;
        return *this;
    }

    /** Logical */
    /**
     * @brief ==
     * @param in: another Interval
     * @return true if both part match
     **/
    bool operator==(const Interval<T> &in)
    {
        return (lower == in.lower) && (upper == in.upper);
    }

    /**
     * @brief !=
     * @param in: another Interval
     * @return true if either part mismatch
     **/
    bool operator!=(const Interval<T> &in)
    {
        return (lower != in.lower) || (upper != in.upper);
    }

    /**
     * @brief ostream
     * @param os
     * @param in: interval
     **/
    friend std::ostream &operator<<(std::ostream &os, const Interval<T> &in)
    {
        os << '[' << in.lower << ", " << in.upper << ']';
        return os;
    }
};

} // namespace Spatial
} // namespace Roahm

#include <eigen3/Eigen/Core>
namespace Eigen
{
template <> struct NumTraits<Roahm::Spatial::Interval<double>>
{
    typedef Roahm::Spatial::Interval<double> Real;
    typedef Roahm::Spatial::Interval<double> NonInteger;
    typedef Roahm::Spatial::Interval<double> Nested;
    typedef double Literal;
    enum
    {
        IsComplex = 0,
        IsInteger = 0,
        IsSigned = 1,
        RequireInitialization = 1,
        ReadCost = 2,
        AddCost = 6,
        MulCost = 6
    };

    static Real epsilon()
    {
        return Real{std::numeric_limits<double>::epsilon(),
                    std::numeric_limits<double>::epsilon()};
    }

    static Real dummy_precision()
    {
        return Real{1e-12, 1e-12};
    }

    static Real highest()
    {
        return Real{std::numeric_limits<double>::max(),
                    std::numeric_limits<double>::max()};
    }

    static Real lowest()
    {
        return Real{-std::numeric_limits<double>::max(),
                    -std::numeric_limits<double>::max()};
    }

    static int digits10()
    {
        return NumTraits<double>::digits10();
    }
};
} // namespace Eigen
