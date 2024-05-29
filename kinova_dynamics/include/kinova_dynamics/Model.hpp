#pragma once
#include "kinova_dynamics/Spatial.hpp"
#include "kinova_dynamics/spatial/SpatialTransform.hpp"
#include "kinova_utils/Utils.hpp"
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <vector>

namespace Roahm
{
namespace Model
{
using namespace Spatial;
using std::vector;

/**
 * @brief possible joint types for %Joint
 */
enum class JointType : uint8_t
{
    UNKNOWN,
    REVOLUTE,
    CONTINUOUS,
    PRISMATIC,
    FLOATING,
    PLANAR,
    FIXED
}; // namespace Model

/**
 * @brief Link type for spatial model
 * @param name: link name
 * @param parent_id: index of the parent body
 * @param child_ids: indices of the child bodies
 * @param SI: spatial inertia
 * @param T: %Transform from parent frame
 * @param S: %SpatialConstraints specifying rotating/prismatic axis
 * @param type: %JointType
 * @param limit: joint limit, only valid for revolute joint
 * @param damping: damping term linear to velocity
 * @param friction: constantly applied friction
 * @param armature: extra inertia caused by motor armature. For details, see
 *  https://mujoco.readthedocs.io/en/latest/XMLreference.html?highlight=armature
 **/
template <typename U> struct SpatialLink
{
    std::string name;
    size_t parent_id;
    vector<size_t> child_ids;
    SpatialInertia<U> SI;
    /** Joint Attributes */
    SpatialTransform<U> T;   // prev to body transform
    SpatialTransform<U> Xbw; // world to body transform
    SpatialTransform<U> CoM; // center of mass
    SpatialAxis<U> S;        // rotation axis
    JointType type = JointType::FIXED;
    Vector2d limit;
    U damping = 0;
    U friction = 0;
    U armature = 0;

    /** Methods */
  public:
    /**
     * @brief get dynamic transform based on joint input
     * @details the transform is already combined with link transform
     * @param q: joint input of type %U
     **/
    SpatialTransform<U> get_transform(const double &q,
                                      const SpatialAxis<U> &axis) const
    {
        switch (this->type)
        {
        case JointType::FIXED:
            // fixed transform
            return T.inverse();
        default: {
            Mat3<U> R = dir2rot(-U(q), Vec3<U>(axis.ang()));
            Vec3<U> p =
                (Mat3<U>::Identity() - R) * axis.ang().cross(axis.lin());
            p = -R.transpose() * p;
            SpatialTransform<U> xJ =
                SpatialTransform<U>{std::move(R), std::move(p)};
            return xJ * T.inverse();
        }
        }
    }
};

using SLinkd = SpatialLink<double>;

/// forward declaration for use in Impl
template <typename T> struct SpatialModel;

/**
 * @brief stores implementation function of spatial model
 **/
namespace SpatialModelImpl
{
/**
 * @brief parse model for double model
 **/
SpatialModel<double> load_double(const std::string &model_path);

/**
 * @brief parse model and add variance
 **/
SpatialModel<Intd> load_interval(const std::string &model_path,
                                 const double variance);

/**
 * @brief create interval model from normal model
 * @param base_model: normal model
 * @param variance: variance of interval
 **/
SpatialModel<Intd> make_interval_model(const SpatialModel<double> &base_model,
                                       const double variance,
                                       const double fric_var = -1);

/**
 * @brief bias a model by given variance
 * @details new_param = param * (1 + variance)
 * @param model: model to bias
 * @param variance: bias of the model to apply
 **/
void bias_model(SpatialModel<double> &model, const double variance);
} // namespace SpatialModelImpl

/**
 * @brief spatial model type
 **/
template <typename T> struct SpatialModel
{
    std::string name;
    vector<SpatialLink<T>> links;
    uint32_t num_joints;
    uint32_t num_dof;
    Vec3<T> gravity{0, 0, -9.81};
    bool is_featherstone = false;

  public:
    /**
     * @brief attach another child model onto the robot
     * @param child_model: child model to attach
     * @param parent_id: index of the link to be attached to
     * NOTE: will fail if child has joints...
     **/
    void attach(const SpatialModel<T> &child_model, size_t parent_id)
    {
        size_t child_root = this->links.size();
        this->num_dof += child_model.num_dof;
        this->num_joints += child_model.num_joints;
        // append link to vector
        this->links.insert(this->links.end(), child_model.links.begin(),
                           child_model.links.end());
        // modify parent ids and child ids
        links[parent_id].child_ids.push_back(child_root);
        links[child_root].parent_id = parent_id;
        for (size_t i = child_root; i < links.size(); i++)
        {
            if (i != child_root)
                links[i].parent_id += child_root;
            // change child id
            for (auto &child_id : links[i].child_ids)
                child_id += child_root;
            // update base world transform
            if (i == child_root)
            {
                links[i].T = links[links[i].parent_id].CoM * links[i].T;
            }
            links[i].Xbw = links[links[i].parent_id].Xbw * links[i].T;
        }
    }

    void attach(const SpatialModel<T> &child_model,
                const std::string &parent_link_name)
    {
        // find link based on name and attach
        for (std::size_t i = 0; i < links.size(); i++)
        {
            if (links[i].name == parent_link_name)
            {
                attach(child_model, i);
                return;
            }
        }
    }

    /**
     * @brief print info about current model
     **/
    std::ostream &print(bool verbose = false,
                        std::ostream &os = std::cout) const
    {
        using std::endl;
        os << "Number of joints: " << num_joints << endl;
        os << "Number of dof: " << num_dof << endl;
        os << endl;
        os << name << endl;
        print_node(verbose, os, 0, "└──");
        return os;
    }

    friend std::ostream &operator<<(std::ostream &os,
                                    const SpatialModel<T> &model)
    {
        model.print(false, os);
        return os;
    }

    /** Double Model Specific */
    /**
     * @brief load normal robot model
     **/
    template <typename U = T,
              std::enable_if_t<std::is_same_v<U, double>, bool> = true>
    static inline SpatialModel<double> load(const std::string &model_path)
    {
        return SpatialModelImpl::load_double(model_path);
    }

    /**
     * @brief bias current model
     * @param variance: variance to bias
     **/
    template <typename U = T,
              std::enable_if_t<std::is_same_v<U, double>, bool> = true>
    inline void bias(const double variance)
    {
        SpatialModelImpl::bias_model(*this, variance);
    }

    /** Interval Model Specifics */
    /**
     * @brief loader for interval robot model
     **/
    template <typename U = T,
              std::enable_if_t<std::is_same_v<U, Intd>, bool> = true>
    static inline SpatialModel<Intd> load(const std::string &model_path,
                                          const double variance)
    {
        return SpatialModelImpl::load_interval(model_path, variance);
    }

    template <typename U = T,
              std::enable_if_t<std::is_same_v<U, Intd>, bool> = true>
    static inline SpatialModel<Intd> from_double(
        const SpatialModel<double> &model, const double variance,
        const double fric_var = -1)
    {
        return SpatialModelImpl::make_interval_model(model, variance, fric_var);
    }

  protected:
    /**
     * @brief recursive helper function to print each node of the model
     **/
    void print_node(bool verbose, std::ostream &os, size_t idx,
                    const std::string padding, bool is_ending = true) const
    {
        using std::endl;
        auto &link = links[idx];
        static const std::vector<std::string> joint_names{
            "Unknown",  "Revolute", "Continuous", "Prismatic",
            "Floating", "Planar",   "Fixed"};
        os << padding << link.name << ": [" << joint_names[uint(link.type)]
           << "]" << endl;
        std::string pad = padding.substr(0, padding.size() - 9);
        std::string apad = "   ";
        if (!is_ending)
            apad = "│  ";
        if (verbose)
        {
            std::string epad = "│   ";
            if (link.child_ids.empty())
                epad = "    ";
            std::string tpad = pad + apad + epad;
            os << tpad << "Mass: " << link.SI.mass() << endl;
            os << tpad << "Center of mass: " << link.SI.com().transpose()
               << endl;
            os << tpad << "Joint axis: " << link.S << endl;
            os << tpad << "Friction: " << link.friction << endl;
            os << tpad << "Damping: " << link.damping << endl;
            os << tpad << "Armature: " << link.armature << endl;
            auto I_c = link.SI.compute_I_c();
            auto I = link.SI.I();
            os << tpad << "Inertia w.r.t. com: " << endl;
            os << tpad << "   " << I_c.row(0) << endl;
            os << tpad << "   " << I_c.row(1) << endl;
            os << tpad << "   " << I_c.row(2) << endl;
            os << tpad << "Inertia w.r.t. origin: " << endl;
            os << tpad << "   " << I.row(0) << endl;
            os << tpad << "   " << I.row(1) << endl;
            os << tpad << "   " << I.row(2) << endl;
            // print transform
            if (idx != 0)
            {
                os << tpad << "Axis: " << link.S << std::endl;
                os << tpad << "Transform: " << endl;
                os << link.T << endl;
            }
        }
        if (link.child_ids.empty())
            return;
        for (uint i = 0; i < link.child_ids.size() - 1; i++)
        {
            print_node(verbose, os, link.child_ids[i], pad + apad + "├──",
                       false);
        }
        print_node(verbose, os, link.child_ids.back(), pad + apad + "└──",
                   true);
    }
};

using SpatialModeld = SpatialModel<double>;
using SpatialModelInt = SpatialModel<Intd>;
} // namespace Model
} // namespace Roahm
