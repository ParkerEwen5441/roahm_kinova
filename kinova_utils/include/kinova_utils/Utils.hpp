#pragma once
#include <eigen3/Eigen/Eigen>
#include <mutex>
#include <vector>

namespace Roahm
{
namespace Utils
{
/**
 * @brief construct Eigen diagonal matrix fromm std::vector
 * @param diagonal: vector of elements on diagonal
 **/
Eigen::MatrixXd DiagonalMatrixXd(std::vector<double> diagonal);

/**
 * @brief convert SRC type to eigen matrix
 * @details need both SRC type to have data() method
 **/
template <typename SRC, typename DST>
DST stl2Vec(const SRC &src, const size_t size = 0)
{
    size_t l = size;
    if (l == 0)
        l = src.size();
    DST dst;
    if (l != dst.size())
    {
        dst.resize(l);
    }
    std::copy(src.data(), src.data() + l, dst.data());
    return dst;
}

/**
 * @brief disable copy constructor and copy assignment operator for derived
 *class
 * @details usage: inherit from this class and declare it as private
 * @ref
 *https://www.boost.org/doc/libs/1_73_0/libs/core/doc/html/core/noncopyable.html
 **/
class noncopyable
{
  public:
    noncopyable() = default;
    ~noncopyable() = default;
    noncopyable(const noncopyable &) = delete;
    noncopyable &operator=(const noncopyable &) = delete;
};

/**
 * @brief stop signal used in thread
 **/
class thread_stop : public std::exception
{
  public:
    virtual const char *what() const noexcept override
    {
        return "Thread signaled to stop";
    }
};

namespace Signal
{
enum class Signal
{
    Interrupt = 2,
};

class SignalException : public std::runtime_error
{
  public:
    template <Signal SIG> static SignalException make()
    {
        return SignalException(SIG);
    }

    SignalException(Signal signal)
        : std::runtime_error("Signal Exception"), signal(signal)
    {
    }

    /**
     * @brief get signal type
     **/
    inline Signal what_signal()
    {
        return signal;
    }

  private:
    Signal signal;
};

[[noreturn]] inline void InterruptSignalHandler([[maybe_unused]] int _ignored)
{
    throw SignalException::make<Signal::Interrupt>();
}

/**
 * @brief Set up the signal handlers for certain signal
 */
void setupSignalHandlers(Signal sig = Signal::Interrupt,
                         void (*handler)(int) = InterruptSignalHandler);
} // namespace Signal
} // namespace Utils
} // namespace Roahm
