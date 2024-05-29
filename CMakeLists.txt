cmake_minimum_required(VERSION 3.10)
project(roahm_kinova)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Include directories
include_directories(
    kinova_utils/include
    kinova_system/include
    kinova_filter/include
    kinova_dynamics/include
)

# Subdirectories
add_subdirectory(kinova_utils)
add_subdirectory(kinova_system)
add_subdirectory(kinova_filter)
add_subdirectory(kinova_dynamics)