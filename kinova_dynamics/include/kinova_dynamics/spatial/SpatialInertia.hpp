/**
 * @brief Spatial Inertia
 * @file
 **/
#pragma once
#include "EigenTypes.hpp"
#include "Helper.hpp"
#include "Twist.hpp"
#include "kinova_dynamics/spatial/SpatialTransform.hpp"

namespace Roahm
{
namespace Spatial
{
template <typename T> class SpatialInertia
{
  protected:
    Mat3<T> ch;
    Mat3<T> I_bar;
    T m;

  public:
    /**
     * @brief default constructor
     **/
    SpatialInertia() : ch(Mat3<T>::Zero()), I_bar(Mat3<T>::Zero()), m(T(0))
    {
    }

    /**
     * @brief copy ctor
     **/
    SpatialInertia(const SpatialInertia<T> &si)
        : ch(si.ch), I_bar(si.I_bar), m(si.m)
    {
    }

    /**
     * @brief assignment operator
     **/
    SpatialInertia<T> &operator=(const SpatialInertia<T> &si)
    {
        ch = si.ch;
        I_bar = si.I_bar;
        m = si.m;
        return *this;
    }

    /**
     * @brief constructor overload
     * @param I_c: tensor of inertia w.r.t. CoM
     * @param CoM: Center of Mass
     * @param m: mass
     **/
    SpatialInertia(const Mat3<T> &I_c, const Vec3<T> &CoM, T m)
        : ch(m * skew_sym_mat(CoM)), I_bar(I_c - ch * ch / m), m(m)
    {
    }

    /**
     * @brief constructor overload for urdf parser
     * @param ixx:
     * @param ixy:
     * @param ixz:
     * @param iyy:
     * @param iyz:
     * @param izz:
     * @param CoM: Center of Mass
     * @param m: mass
     **/
    SpatialInertia(T ixx, T ixy, T ixz, T iyy, T iyz, T izz, const Vec3<T> &CoM,
                   T m)
        : ch(m * skew_sym_mat(CoM)), I_bar(), m(m)
    {
        Mat3<T> I_c;
        I_c << ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz;
        I_bar = I_c - ch * ch / m;
    }

    /**
     * @breif get the inertia matrix w.r.t. origin of the frame
     * @return tensor of inertia
     **/
    const Mat3<T> &I() const
    {
        return I_bar;
    }

    /**
     * @brief get the inertia matrix w.r.t. CoM
     * @details will be slow to call this since it's not directly stored
     * @return tensor of inertia matrix w.r.t. CoM
     **/
    Mat3<T> compute_I_c() const
    {
        return I() + m * ch * ch;
    }

    /**
     * @brief get mass
     * @return mass value
     **/
    inline T mass() const
    {
        return m;
    }

    /**
     * @brief get center of mass
     * @return const reference of center of mass vector
     **/
    Vec3<T> com() const
    {
        Mat3<T> ch_n = ch / mass();
        Vec3<T> com(-ch_n(1, 2), ch_n(0, 2), -ch_n(0, 1));
        return com;
    }

    const Mat3<T> &m_c_hat() const
    {
        return ch;
    }

    /**
     * @brief apply inertia on twist
     * @param v: %Twist to apply
     * @return general force of type %Wrench
     **/
    Wrench<T> operator*(const Twist<T> &v) const
    {
        return Wrench<T>{I() * v.angular + ch * v.linear,
                         mass() * v.linear - ch * v.angular};
    }

    /**
     * @brief add two %SpatialInertia together
     **/
    SpatialInertia<T> &operator+=(const SpatialInertia<T> &si)
    {
        I_bar += si.I_bar;
        ch += si.ch;
        m += si.m;
        return *this;
    }
    SpatialInertia<T> operator+(const SpatialInertia<T> &si)
    {
        SpatialInertia<T> res(*this);
        res += si;
        return res;
    }

    /**
     * @brief apply transform on spatial inertia
     **/
    SpatialInertia<T> apply_transform(const SpatialTransform<T> &X)
    {
        Mat3<T> p_hat = skew_sym_mat(X.translation());

        const Mat3<T> &R = X.rotation();
        Mat3<T> mRp_hat = m * R * p_hat;
        SpatialInertia<T> new_si;
        new_si.ch = R * ch * R.transpose() - m * R * p_hat * R.transpose();
        new_si.I_bar =
            (R * (I_bar + 2 * ch * p_hat) - mRp_hat * p_hat) * R.transpose();
        new_si.m = m;
        return new_si;
    }

  protected:
};
} // namespace Spatial
} // namespace Roahm
