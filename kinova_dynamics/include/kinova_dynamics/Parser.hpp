#pragma once
#include "kinova_dynamics/Model.hpp"
#include <string>
#include <vector>

namespace Roahm
{
namespace Model
{
class Parser
{
  public:
    /**
     * @brief parse model
     * @param model_path: path to model description file
     * @param model_format: format of the model, if not specified, will
     *        determine according to file extension.
     *        Options: urdf, mjcf
     **/
    static SpatialModel<double> parse(const std::string &model_path,
                                      std::string model_format = "");

  protected:
    /**
     * @brief parse model from urdf file
     * @param urdf_path: path to urdf file
     **/
    static SpatialModel<double> parse_urdf(const std::string &urdf_path);

    /**
     * @brief parse model from mjcf file
     * @param mjcf_path: path to mjcf file
     **/
    static SpatialModel<double> parse_mjcf(const std::string &mjcf_path);
};
} // namespace Model
} // namespace Roahm
