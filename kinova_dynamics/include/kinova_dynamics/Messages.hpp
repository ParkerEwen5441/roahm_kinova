/**
 * @brief Message structs and classes
 * @file
 **/
#pragma once
#include <chrono>
#include <cstdlib>
#include <eigen3/Eigen/Eigen>

namespace Roahm
{
namespace msgs
{
/**
 * @warning default Ctor is provided because Port need them to initialize
 *          should explicitly set length if using elsewhere
 **/
struct Measurement
{
    std::size_t frame_id;
    Eigen::VectorXd pos;
    Eigen::VectorXd vel;
    Eigen::VectorXd torque;
    decltype(std::chrono::system_clock::now()) timestamp;
    /** Ctors */
    Measurement() = default;
    Measurement(uint32_t len)
        : frame_id(0), pos(len), vel(len), torque(len),
          timestamp(std::chrono::system_clock::now())
    {
        pos.setZero();
        vel.setZero();
        torque.setZero();
    }
};

struct ControlMsg
{
    std::size_t frame_id;
    Eigen::VectorXd control;
    bool is_gripper_open;
    /** Ctors */
    ControlMsg() = default;
    ControlMsg(uint32_t len) : frame_id(0), control(len), is_gripper_open(true)
    {
        control.setZero();
    }
};

/**
 * @brief trajectory to be passed to controllerblock
 **/
struct AccTrajectory
{
    Eigen::VectorXd k;
    Eigen::VectorXd acc;
    Eigen::VectorXd vel;
    Eigen::VectorXd pos;
    double start_time;
    double duration;
    double trajectory_duration;
    bool is_gripper_open;
    bool is_bernstein;
    bool reset;

    /** Ctors */
    AccTrajectory() = default;
    AccTrajectory(uint32_t len)
        : k(len), acc(len), vel(len), pos(len), start_time(0), duration(0.5), trajectory_duration(1.0),
          is_gripper_open(false), is_bernstein(false), reset(false)
    {
        k.setZero();
        acc.setZero();
        vel.setZero();
        pos.setZero();
    }
};
} // namespace msgs
} // namespace Roahm
