#include "kinova_dynamics/Parser.hpp"
#include "kinova_dynamics/Model.hpp"
#include "kinova_dynamics/Spatial.hpp"
#include "kinova_dynamics/spatial/SpatialAxis.hpp"
#include "kinova_utils/Logger.hpp"
#include "kinova_utils/Utils.hpp"
#include <bits/stdint-uintn.h>
#include <cassert>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <sstream>
#include <stack>
#include <string>
#include <tinyxml2.h>
#include <vector>

Roahm::Logger logger("parser");
/** HELPER FUNCTIONS */
/**
 * @brief output the given error and exit
 **/
void xml_error(const std::string &msg)
{
    logger.error(msg);
    exit(-1);
}

/**
 * @brief convert standard model read from description files to roy featherstone
 *style
 * @param rob: robot model to convert to
 *
 **/
void featherstone_conversion(Roahm::Model::SpatialModel<double> &rob)
{
    if (rob.is_featherstone)
        return;

    Roahm::Model::STransformd Xwj = rob.links[0].T;
    for (uint32_t i = 0; i < rob.links.size(); i++)
    {
        auto &link = rob.links[i];
        // Convert S from joint frame to world frame
        const std::size_t pid = link.parent_id;

        if (pid != std::string::npos)
        {
            Xwj = link.T * Xwj;
        }

        Roahm::Model::Twistd s = Xwj.inv_apply(link.S.as_twist());
        link.S = Roahm::Spatial::SpatialAxis<double>(s.angular, s.linear);

        // Convert I from joint frame to body CoM frame
        link.SI = link.SI.apply_transform(link.CoM);

        // Convert Xtree from joint-to-joint to CoM-to-CoM (backwards)
        Roahm::Model::STransformd prevCoM;
        if (pid != std::string::npos)
        {
            prevCoM = rob.links[pid].CoM;
        }

        link.T = prevCoM * link.T.inverse() * link.CoM.inverse();

        if (pid != std::string::npos)
        {
            link.Xbw = rob.links[pid].Xbw * link.T;
        }
        else
        {
            link.Xbw = link.T;
        }
    }
    rob.is_featherstone = true;
}

namespace xml
{
using namespace tinyxml2;
using std::string;

/**
 * @brief check whether the attribute exist
 * @param elm: element to check
 * @param attr: attribute name
 * @return if the attribute exists
 **/
bool has_attribute(const XMLElement *elm, const string &attr)
{
    return (elm->Attribute(attr.c_str()) != nullptr);
}

/**
 * @brief get attribute of element
 * @param elm: element to check
 * @param attr: attribute name
 * @return empty string if no such attribute, else just return value
 **/
string get_attribute(const XMLElement *elm, const string &attr)
{
    const char *cres = elm->Attribute(attr.c_str());
    if (cres == nullptr)
        return "";
    return string(cres);
}

/**
 * @brief parse attributes as array of double
 * @return vector of double storing parsed result
 **/
std::vector<double> parse_array_attr(const XMLElement *elm,
                                     const std::string &attr,
                                     const size_t size = 0,
                                     const char separator = ' ')
{
    // get attribute value in str
    std::string val = get_attribute(elm, attr);
    if (val.empty())
    {
        xml_error("Empty attribute");
    }
    std::stringstream ss(val);
    string str;
    std::vector<double> res;
    // separate string using separator
    while (std::getline(ss, str, separator))
    {
        // convert to double
        res.push_back(std::stod(str));
    }
    // check size
    if (size != 0)
        assert(res.size() == size);
    return res;
}
} // namespace xml

namespace mjcf
{
const auto toVec3 = [](const std::vector<double> &vec) {
    return Roahm::Utils::stl2Vec<std::vector<double>, Eigen::Vector3d>(vec, 3);
};

/**
 * @return either attr value or default if there's no
 **/
std::string attr_or_default(tinyxml2::XMLElement *elm, const std::string &attr,
                            const std::string &def_val)
{
    if (!xml::has_attribute(elm, attr))
        return def_val;
    return xml::get_attribute(elm, attr);
}

/**
 * @brief get rotation matrix from xml element
 * @param body: xml element to search
 * @param required: if required and not provided, trigger error
 * @return rotation matrix or Identity if not provided
 **/
Eigen::Matrix3d get_rotation(tinyxml2::XMLElement *body, bool required = false)
{
    using namespace Eigen;
    Matrix3d rot = Matrix3d::Identity();
    // check types of rotation 1 by 1
    if (xml::has_attribute(body, "quat"))
    {
        auto qvec = xml::parse_array_attr(body, "quat");
        rot =
            Quaterniond(qvec[0], qvec[1], qvec[2], qvec[3]).toRotationMatrix();
    }
    else if (xml::has_attribute(body, "euler"))
    {
        auto evec = xml::parse_array_attr(body, "euler");
        rot = (AngleAxisd(evec[0], Vector3d::UnitX()) *
               AngleAxisd(evec[1], Vector3d::UnitY()) *
               AngleAxisd(evec[2], Vector3d::UnitZ()))
                  .toRotationMatrix();
    }
    else if (required)
    {
        xml_error("Rotation for " + std::string(body->Name()) +
                  " required but not given");
    }
    return rot;
}

/**
 * @brief get spatial transform matrix from body
 * @param body: xml body to get data
 * @return %SpatialTransform
 **/
Roahm::Spatial::STransformd get_transform(tinyxml2::XMLElement *body)
{
    using namespace Eigen;
    using Roahm::Spatial::STransformd;
    // translation
    Vector3d p = Vector3d::Zero();
    if (xml::has_attribute(body, "pos"))
        p = toVec3(xml::parse_array_attr(body, "pos"));
    // rotation and return
    return STransformd{get_rotation(body).transpose(), std::move(p)};
}

/**
 * @brief get inertia from mjcf field
 * @param inertia: inertia xml tag, assumed to be valid
 * @return %SpatialInertia
 **/
Roahm::Spatial::SpatialInertia<double> get_inertia(
    tinyxml2::XMLElement *inertial)
{
    using namespace Eigen;
    // center of mass
    Vector3d CoM = toVec3(xml::parse_array_attr(inertial, "pos"));
    // mass
    double mass = std::stod(xml::get_attribute(inertial, "mass"));
    // moment of inertia
    std::vector<double> i;
    if (xml::has_attribute(inertial, "diaginertia"))
    {
        i = xml::parse_array_attr(inertial, "diaginertia");
        i.resize(6, 0);
    }
    else if (xml::has_attribute(inertial, "fullinertia"))
    {
        i = xml::parse_array_attr(inertial, "fullinertia");
    }
    else
    {
        xml_error("inertia should have either diaginertia or "
                  "fullinertia field!");
    }
    Matrix3d I_c;
    I_c << i[0], i[3], i[4], i[3], i[1], i[5], i[4], i[5], i[2];
    // get rotation
    Matrix3d rot = get_rotation(inertial);
    return Roahm::Spatial::SpatialInertia<double>(rot * I_c * rot.transpose(),
                                                  CoM, mass);
}
} // namespace mjcf

/** ROAHM */
namespace Roahm
{
namespace Model
{
SpatialModel<double> Parser::parse(const std::string &model_path,
                                   std::string model_format)
{
    // check model format automatically
    if (model_format.empty())
    {
        size_t dot = model_path.find_last_of('.');
        model_format = model_path.substr(dot + 1);
    }

    SpatialModel<double> rob;
    // parse depends on type
    if (model_format == "mjcf")
    {
        // logger.info("Parsing Mujoco Description file...");
        rob = Parser::parse_mjcf(model_path);
    }
    else if (model_format == "urdf")
    {
        rob = Parser::parse_urdf(model_path);
    }
    else
    {
        xml_error("model format [" + model_format + "] not supported");
        exit(-1);
    }

    // Convert to Roy Featherstone style model
    featherstone_conversion(rob);

    return rob;
}

SpatialModel<double> Parser::parse_mjcf(const std::string &mjcf_path)
{
    using namespace mjcf;
    using namespace tinyxml2;

    // load file
    XMLDocument mjcf;
    mjcf.LoadFile(mjcf_path.c_str());
    if (mjcf.Error())
    {
        logger.error("Unable to parse mujoco description file \"{}\"",
                     mjcf_path);
        mjcf.PrintError();
        exit(1);
    }

    // get root element
    XMLElement *root = mjcf.FirstChildElement("mujoco")
                           ->FirstChildElement("worldbody")
                           ->FirstChildElement("body");
    // ensure no sibling
    {
        XMLElement *bd_sibling = root->NextSiblingElement("body");
        if (bd_sibling != nullptr)
        {
            xml_error("parser only support description file with "
                      "only 1 root body element!");
        }
    }
    SpatialModel<double> rob;
    rob.name = attr_or_default(mjcf.FirstChildElement("mujoco"), "model",
                               "unknown robot");

    // create stack for depth first search
    size_t npos = std::string::npos;
    std::stack<std::pair<XMLElement *, size_t>> link_stack;
    link_stack.push({root, npos});

    /** LOOP */
    while (!link_stack.empty())
    {
        // get link to process
        auto link_info = link_stack.top();
        link_stack.pop();
        XMLElement *body = link_info.first;
        const size_t pid = link_info.second; // parent id
        // handle link
        SpatialLink<double> link;
        link.name = attr_or_default(body, "name", "unknown link");
        // setting id and also add for parent
        link.parent_id = pid;
        size_t id = rob.links.size();
        /** Inertia */
        XMLElement *inertial = body->FirstChildElement("inertial");
        // skip virtual link
        if (inertial == nullptr ||
            attr_or_default(inertial, "mass", "0") == "0")
            continue;
        link.SI = get_inertia(inertial);
        // center of mass
        link.CoM = STransformd::pure_translation(
            toVec3(xml::parse_array_attr(inertial, "pos")));

        /** Joint */
        // transform
        link.T = get_transform(body);
        if (pid != npos)
        {
            // register as child
            rob.links[pid].child_ids.push_back(id);

            // joint informations
            XMLElement *joint = body->FirstChildElement("joint");
            if (joint != nullptr)
            {
                bool limited = (xml::get_attribute(joint, "limited") == "true");
                std::string stype;

                // joint type (default: hinge)
                stype = attr_or_default(joint, "type", "hinge");

                if (stype == "hinge")
                {
                    // TODO: consider switch to mujoco type which is more
                    // consistent
                    if (limited)
                        link.type = JointType::REVOLUTE;
                    else
                        link.type = JointType::CONTINUOUS;
                }
                else if (stype == "slide")
                {
                    link.type = JointType::PRISMATIC;
                }
                else
                {
                    // TODO adding ball joint
                    logger.warn("Link type {} not supported yet", stype);
                    link.type = JointType::FIXED;
                }

                // angle/movement limit
                if (limited)
                {
                    auto vrange = xml::parse_array_attr(joint, "range");
                    link.limit = Utils::stl2Vec<std::vector<double>, Vector2d>(
                        vrange, 2);
                }

                // joint axis
                if (!xml::has_attribute(joint, "axis"))
                    xml_error("Joint no axis!");
                // TODO add ball joint
                if (stype == "hinge" || stype == "slide")
                {
                    link.S =
                        SAxisd(toVec3(xml::parse_array_attr(joint, "axis")),
                               stype[0] == 's');
                }

                // friction
                link.friction =
                    std::stod(attr_or_default(joint, "frictionloss", "0"));

                // damping
                link.damping =
                    std::stod(attr_or_default(joint, "damping", "0"));

                // armature
                link.armature =
                    std::stod(attr_or_default(joint, "armature", "0"));
            }
        }
        // push link
        rob.links.push_back(std::move(link));
        XMLElement *child = body->FirstChildElement("body");
        while (child != nullptr)
        {
            link_stack.push({child, id});
            child = child->NextSiblingElement("body");
        }
    }
    rob.num_joints = rob.links.size() - 1;
    rob.num_dof = 0;
    for (auto &link : rob.links)
    {
        if (link.type != JointType::FIXED)
        {
            rob.num_dof++;
        }
    }

    return rob;
}

} // namespace Model
} // namespace Roahm

/** URDF Parser */
/// TODO: Change to tinyxml to reduce dependency
#ifdef URDF_SUPPORT
#include <urdf/model.h>
#include <urdf_model/types.h>
namespace urdf
{
/**
 * @brief get transform from urdf::Pose
 * @param pose: %Pose type
 **/
Roahm::Spatial::STransformd getTransform(const urdf::Pose &pose)
{
    using namespace Eigen;
    using namespace Roahm::Spatial;
    STransformd tr(Quaterniond(pose.rotation.w, pose.rotation.x,
                               pose.rotation.y, pose.rotation.z)
                       .toRotationMatrix()
                       .transpose(),
                   Vector3d(pose.position.x, pose.position.y, pose.position.z));
    return tr;
}
} // namespace urdf

Roahm::Model::SpatialModel<double> Roahm::Model::Parser::parse_urdf(
    const std::string &urdf_path)
{
    using namespace urdf;
    static const size_t npos = std::string::npos;
    /** PREPARE */
    urdf::Model model;
    if (!model.initFile(urdf_path))
    {
        xml_error("Failed to parse urdf file");
    }
    SpatialModel<double> rob;
    rob.name = model.getName();
    // create stack for DFS
    std::stack<std::pair<urdf::LinkSharedPtr, size_t>> link_stack;
    link_stack.push({model.root_link_, npos});
    // base transformation to be applied
    STransformd base_tran;

    /** LOOP */
    while (!link_stack.empty())
    {
        // get link to process
        auto link_info = link_stack.top();
        link_stack.pop();
        const auto &urdf_link = link_info.first;
        const size_t pid = link_info.second; // parent id

        /** check virtual link */
        if (urdf_link->inertial == nullptr)
        {
            // FIXME: root link can't be virtual I guess...
            if (urdf_link->parent_joint == nullptr)
            {
                xml_error("root link no inertial");
            }
            // virtual and last, clear base_tran
            if (urdf_link->child_joints.empty())
            {
                base_tran.setIdentity();
                continue;
            }
            base_tran =
                getTransform(
                    urdf_link->parent_joint->parent_to_joint_origin_transform) *
                base_tran;
            // pushing children
            for (auto joint : urdf_link->child_joints)
            {
                urdf::LinkSharedPtr next;
                model.getLink(joint->child_link_name, next);
                link_stack.push({next, pid});
            }
            continue;
        }
        /** must be real link now */
        SpatialLink<double> link;
        link.name = urdf_link->name;
        // setting id and also add for parent
        link.parent_id = pid;
        size_t id = rob.links.size();
        // parent joints
        if (urdf_link->parent_joint != nullptr)
        {
            // pid
            rob.links[pid].child_ids.push_back(id);
            auto joint = urdf_link->parent_joint;
            link.T = getTransform(joint->parent_to_joint_origin_transform);
            link.type = (JointType)joint->type;
            switch (link.type)
            {
            case JointType::REVOLUTE:
                // limit for revolute joint
                link.limit =
                    Vector2d(joint->limits->lower, joint->limits->upper);
                [[fallthrough]];
            case JointType::CONTINUOUS:
                link.S = SAxisd(
                    Vector3d(joint->axis.x, joint->axis.y, joint->axis.z),
                    false);
                break;
            case JointType::PRISMATIC:
                link.S = SAxisd(
                    Vector3d(joint->axis.x, joint->axis.y, joint->axis.z),
                    true);
                break;
            case JointType::FIXED:
                break;
            default:
                xml_error("Joint type not supported yet");
            }
        }
        auto i = urdf_link->inertial;
        Vector3d CoM(i->origin.position.x, i->origin.position.y,
                     i->origin.position.z);
        link.SI = SpatialInertia<double>(i->ixx, i->ixy, i->ixz, i->iyy, i->iyz,
                                         i->izz, CoM, i->mass);
        link.CoM = STransformd::pure_translation(std::move(CoM));

        // add link
        rob.links.push_back(std::move(link));
        for (auto joint : urdf_link->child_joints)
        {
            urdf::LinkSharedPtr next;
            model.getLink(joint->child_link_name, next);
            link_stack.push({next, id});
        }
    }
    rob.num_joints = rob.links.size() - 1;
    rob.num_dof = 0;
    for (auto &link : rob.links)
    {
        if (link.type != JointType::FIXED)
        {
            rob.num_dof++;
        }
    }
    return rob;
}
// #else
// /// throw not implemented error
// Roahm::Model::SpatialModel<double> Roahm::Model::Parser::parse_urdf([[maybe_unused]] const std::string &urdf_path)
// {
//     throw std::runtime_error("URDF Support not built");
// }
// #endif
