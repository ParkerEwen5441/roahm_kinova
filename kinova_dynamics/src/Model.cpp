#include "kinova_dynamics/Model.hpp"
#include "kinova_dynamics/Parser.hpp"
#include "kinova_dynamics/Spatial.hpp"
#include <eigen3/Eigen/src/Core/util/XprHelper.h>
#include <stdexcept>

namespace Roahm
{
namespace Model
{
/**
 * @brief variate base value by variance
 * @return interval with specified variation
 **/
Intd variate(const double base, const double variance)
{
    return Intd((1 - variance) * base, (1 + variance) * base);
}

namespace SpatialModelImpl
{
SpatialModel<double> load_double(const std::string &model_path)
{
    return Parser::parse(model_path);
}

SpatialModel<Intd> load_interval(const std::string &model_path,
                                 const double variance)
{
    auto base_model = Parser::parse(model_path);
    return make_interval_model(base_model, variance);
}

/// TODO: Neeeeeeeeed to verify more than this
SpatialModel<Intd> make_interval_model(const SpatialModel<double> &base_model,
                                       const double variance,
                                       const double fric_var)
{
    // check variance input
    if (variance < 0)
        throw std::runtime_error(
            "interval model variance should be non-negative!");

    auto num_link = base_model.links.size();
    SpatialModel<Intd> model;
    /// Basic infos
    model.gravity = Vec3<Intd>(base_model.gravity(0), base_model.gravity(1),
                               base_model.gravity(2));
    model.name = base_model.name;
    model.num_dof = base_model.num_dof;
    model.num_joints = base_model.num_joints;
    /// links
    for (std::size_t i = 0; i < num_link; i++)
    {
        const auto &base_link = base_model.links[i];
        SpatialLink<Intd> link;
        // non-variable part
        link.limit = base_link.limit;
        link.type = base_link.type;
        link.name = base_link.name;
        link.parent_id = base_link.parent_id;
        link.child_ids = base_link.child_ids;
        // variable part
        // inertia
        {
            auto I_c = base_link.SI.compute_I_c();
            auto com = base_link.SI.com();
            Intd m_int = variate(base_link.SI.mass(), variance);
            Mat3<Intd> I_c_int;
            for (uint i = 0; i < 9; i++)
                I_c_int(i) = variate(I_c(i), variance);
            Vec3<Intd> com_int(variate(com[0], variance),
                               variate(com[1], variance),
                               variate(com[2], variance));
            link.SI = IntSInertiad(I_c_int, com_int, m_int);
        }
        // transform and axis
        {
            auto s_ang = Vec3<Intd>(variate(base_link.S.ang()[0], 0.),
                                    variate(base_link.S.ang()[1], 0.),
                                    variate(base_link.S.ang()[2], 0.));

            auto s_lin = Vec3<Intd>(variate(base_link.S.lin()[0], 0.),
                                    variate(base_link.S.lin()[1], 0.),
                                    variate(base_link.S.lin()[2], 0.));
            link.S = IntSAxisd(s_ang, s_lin);
        }
        {
            auto &R = base_link.T.rotation();
            auto &t = base_link.T.translation();
            Vec3<Intd> t_int(variate(t(0), 0.), variate(t(1), 0.),
                             variate(t(2), 0.));
            Mat3<Intd> R_int;
            for (uint i = 0; i < 9; i++)
                R_int(i) = variate(R(i), 0.);
            link.T = IntSTransformd(std::move(R_int), std::move(t_int));
        }
        {
            auto &R = base_link.Xbw.rotation();
            auto &t = base_link.Xbw.translation();
            Vec3<Intd> t_int(variate(t(0), 0.), variate(t(1), 0.),
                             variate(t(2), 0.));
            Mat3<Intd> R_int;
            for (uint i = 0; i < 9; i++)
                R_int(i) = variate(R(i), 0.);
            link.Xbw = IntSTransformd(std::move(R_int), std::move(t_int));
        }
        {
            auto &R = base_link.CoM.rotation();
            auto &t = base_link.CoM.translation();
            Vec3<Intd> t_int(variate(t(0), 0.), variate(t(1), 0.),
                             variate(t(2), 0.));
            Mat3<Intd> R_int;
            for (uint i = 0; i < 9; i++)
                R_int(i) = variate(R(i), 0.);
            link.CoM = IntSTransformd(std::move(R_int), std::move(t_int));
        }

        // friction, damping, armature
        if (fric_var != -1)
            link.friction = variate(base_link.friction, fric_var);
        else
            link.friction = variate(base_link.friction, variance);

        link.damping = variate(base_link.damping, variance);
        link.armature = variate(base_link.armature, variance);
        model.links.push_back(link);
    }
    return model;
}

void bias_model(SpatialModel<double> &model, const double variance)
{
    const auto bias = [](const double base, const double variance) {
        return base * (1 + variance);
    };
    for (auto &link : model.links)
    {
        // inertia
        Mat3<double> I_c = link.SI.compute_I_c();
        Vec3<double> com = link.SI.com();
        double m = bias(link.SI.mass(), variance);
        for (uint i = 0; i < 9; i++)
            I_c(i) = bias(I_c(i), variance);
        com << bias(com[0], variance), bias(com[1], variance),
            bias(com[2], variance);
        link.SI = SInertiad(I_c, com, m);
        // transform and axis
        {
            auto s_ang = Vec3<double>(bias(link.S.ang()[0], variance),
                                      bias(link.S.ang()[1], variance),
                                      bias(link.S.ang()[2], variance));

            auto s_lin = Vec3<double>(bias(link.S.lin()[0], variance),
                                      bias(link.S.lin()[1], variance),
                                      bias(link.S.lin()[2], variance));
            link.S = SAxisd(s_ang, s_lin);
        }
        {
            auto &R = link.T.rotation();
            auto &t = link.T.translation();
            Vec3<double> t_biased(bias(t(0), variance), bias(t(1), variance),
                                  bias(t(2), variance));
            Mat3<double> R_biased;
            for (uint i = 0; i < 9; i++)
                R_biased(i) = bias(R(i), variance);
            link.T = STransformd(std::move(R_biased), std::move(t_biased));
        }
        {
            auto &R = link.Xbw.rotation();
            auto &t = link.Xbw.translation();
            Vec3<double> t_int(bias(t(0), variance), bias(t(1), variance),
                               bias(t(2), variance));
            Mat3<double> R_int;
            for (uint i = 0; i < 9; i++)
                R_int(i) = bias(R(i), variance);
            link.Xbw = STransformd(std::move(R_int), std::move(t_int));
        }
        {
            auto &R = link.CoM.rotation();
            auto &t = link.CoM.translation();
            Vec3<double> t_int(bias(t(0), variance), bias(t(1), variance),
                               bias(t(2), variance));
            Mat3<double> R_int;
            for (uint i = 0; i < 9; i++)
                R_int(i) = bias(R(i), variance);
            link.CoM = STransformd(std::move(R_int), std::move(t_int));
        }

        // friction, damping, armature
        link.friction = bias(link.friction, variance);
        link.damping = bias(link.damping, variance);
        link.armature = bias(link.armature, variance);
    }
}

} // namespace SpatialModelImpl
} // namespace Model
} // namespace Roahm
