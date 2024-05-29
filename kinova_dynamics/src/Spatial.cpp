#include "kinova_dynamics/Spatial.hpp"

namespace Roahm
{
namespace Spatial
{
/**
 * @brief implement math functions
 **/
Interval<double> sin(Interval<double> val)
{
    return Interval<double>{std::sin(val.lower), std::sin(val.upper)};
}

Interval<double> cos(Interval<double> val)
{
    return Interval<double>{std::cos(val.lower), std::cos(val.upper)};
}
} // namespace Spatial
} // namespace Roahm
