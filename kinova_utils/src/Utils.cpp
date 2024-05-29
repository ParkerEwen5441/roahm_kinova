#include "kinova_utils/Utils.hpp"
#include <cstddef>
#include <errno.h>
#include <iostream>
#include <signal.h>
#include <stdexcept>

namespace Roahm
{
namespace Utils
{
Eigen::MatrixXd DiagonalMatrixXd(std::vector<double> diagonal)
{
    using namespace Eigen;
    return Map<VectorXd, 0, InnerStride<1>>(&diagonal[0], diagonal.size())
        .asDiagonal();
}

namespace Signal
{
void setupSignalHandlers(Signal sig, void (*handler)(int))
{
    if (signal((int)sig, handler) == SIG_ERR)
    {
        throw std::runtime_error("Error setting up signal handlers");
    }
}
} // namespace Signal
} // namespace Utils
} // namespace Roahm
