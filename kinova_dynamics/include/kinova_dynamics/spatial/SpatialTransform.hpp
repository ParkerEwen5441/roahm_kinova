/**
 * @brief spatial homogeneous transform class
 * @file
 **/
#pragma once
#include "EigenTypes.hpp"
#include "Twist.hpp"
#include "SpatialAxis.hpp"
#include "Wrench.hpp"
#include "kinova_dynamics/spatial/Helper.hpp"
#include <ostream>

namespace Roahm
{
namespace Spatial
{
using namespace Eigen;
template <typename T> class SpatialTransform
{
  protected:
    // rotation and translation part
    Mat3<T> _R;
    Vec3<T> _p;

  public:
    /**
     * @brief create a transformation of pure translation
     * @param translation: translation
     **/
    static SpatialTransform<T> pure_translation(Vec3<T> &&pos)
    {
        return SpatialTransform<T>{Mat3<T>::Identity(), std::move(pos)};
    }

    /**
     * @brief create a transformation of pure rotation
     * @param rot: rotation matrix
     **/
    static SpatialTransform<T> pure_rotation(Mat3<T> &&rot)
    {
        return SpatialTransform<T>{std::move(rot), Vec3<T>::Zero()};
    }

    /**
     * @brief construct from euler angle
     * @param euler: euler angle
     * @param pos: translation
     **/
    static SpatialTransform<T> from_euler(const Vec3<T> &euler,
                                          const Vec3<T> &pos)
    {
        return SpatialTransform<T>{(AngleAxis<T>(euler[0], Vector3d::UnitX()) *
                                    AngleAxis<T>(euler[1], Vector3d::UnitY()) *
                                    AngleAxis<T>(euler[2], Vector3d::UnitZ()))
                                       .toRotationMatrix(),
                                   pos};
    }

    /**
     * @brief default constructor setting transform to identity
     **/
    SpatialTransform() : _R{Mat3<T>::Identity()}, _p{Vec3<T>::Zero()}
    {
    }

    /**
     * @brief constructor from parts
     * @details enforce rvalue to be passed
     * @param R: rvalue ref to R
     * @param p: rvalue ref to p
     **/
    SpatialTransform(Mat3<T> &&R, Vec3<T> &&p)
        : _R{std::move(R)}, _p{std::move(p)}
    {
    }

    // TODO add cons(constraint, angle)

    /**
     * @brief constructor from quaternion and translation
     * @param quat: quaternion number
     * @param pos: %Vector of translation
     **/
    SpatialTransform(const Quaternion<T> &quat, const Vec3<T> &pos)
        : _R{quat.toRotationMatrix()}, _p{pos}
    {
    }

    /** Methods */
    /**
     * @brief get rotation part
     **/
    const Mat3<T> &rotation() const
    {
        return _R;
    }

    /**
     * @brief get translation part
     **/
    const Vec3<T> &translation() const
    {
        return _p;
    }

    /**
     * @brief calculate inverse transform
     * @return inverse of current transformation
     **/
    SpatialTransform<T> inverse() const
    {
        return SpatialTransform<T>{_R.transpose(), -_R * _p};
    }

    /**
     * @brief set translation part
     **/
    void set_p(const Vec3<T> &p)
    {
        _p = p;
    }

    /**
     * @brief set rotation matrix
     **/
    void set_R(const Mat3<T> &R)
    {
        _R = R;
    }

    /**
     * @brief set the transform to identity
     **/
    void setIdentity()
    {
        _R = Mat3<T>::Identity();
        _p = Vec3<T>::Zero();
    }

    /** Operators */
    /**
     * @brief concat transforms to form a new one
     * @param st: anotheer %SpatialTransform
     * @return new transform
     **/
    SpatialTransform<T> operator*(const SpatialTransform<T> &st) const
    {
        return SpatialTransform<T>{_R * st._R, st._p + st._R.transpose() * _p};
    }

    /**
     * @brief apply transform on %Twist
     * @details Effi.pdf of Featherstone's lecture notes
     * @param v: %Twist to be applied
     * @return %Twist with transform applied
     **/
    Twist<T> apply(const Twist<T> &v) const
    {
        /** Featherstone */
        return Twist<T>{_R * v.angular, _R * (v.linear - _p.cross(v.angular))};
        /** Mine */
        // return Twist<T>{_R * v.angular,
        //                 _R * v.linear + _p.cross(_R * v.angular)};
    }

    /**
     * @brief apply inverse transform on %Twist
     * @details Effi.pdf of Featherstone's lecture notes
     * @param v: %Twist to be applied
     * @return %Twist with transform applied
     **/
    Twist<T> inv_apply(const Twist<T> &v) const
    {
        /** Featherstone */
        Vec3<T> new_w = _R.transpose() * v.angular;
        return Twist<T>{new_w, _R.transpose() * v.linear + _p.cross(new_w)};
        /** Mine */
        // return Twist<T>{_R.transpose() * v.angular,
        //                 _R.transpose() * (v.linear - _p.cross(v.angular))};
    }

    /**
     * @brief apply f1 = X*f2 on Wrench, where X* = X^(-T)
     * @details Effi.pdf of Featherstone's lecture notes
     * @param w: %Wrench to be applied
     * @return %Wrench with transform applied
     **/
    Wrench<T> T_apply(const Wrench<T> &w) const
    {
        /** Featherstone */
        return Wrench<T>{_R.transpose() * w.angular +
                             _p.cross(_R.transpose() * w.linear),
                         _R.transpose() * w.linear};
        /** Mine */
        // return Wrench<T>{_R.transpose() * (w.angular - _p.cross(w.linear)),
        //                  _R.transpose() * w.linear};
    }

    /**
     * @brief apply inverse of X* on spatial vector
     * @details Effi.pdf of Featherstone's lecture notes
     * @param w: %Wrench to be applied
     * @return %Wrench with transform applied
     **/
    Wrench<T> inv_T_apply(const Wrench<T> &w) const
    {
        /** Featherstone */
        return Wrench<T>{_R * (w.angular - _p.cross(w.linear)), _R * w.linear};
        /** Mine */
        // return Wrench<T>{_R * w.angular + _p.cross(_R * w.linear),
        //                  _R * w.linear};
    }

    /**
     * @brief output overload
     * @param os: ostream
     * @param t: tansform
     **/
    friend std::ostream &operator<<(std::ostream &os,
                                    const SpatialTransform<T> &t)
    {
        os << "R: " << std::endl << t._R << std::endl;
        os << "p: " << std::endl << t._p.transpose();
        return os;
    }
};

} // namespace Spatial
} // namespace Roahm
