cmake_minimum_required(VERSION 3.10)
project(roahm_kinova)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Find required packages
find_package(Eigen3 REQUIRED)
find_package(TinyXML2 REQUIRED)
find_package(fmt REQUIRED)

# Include directories
include_directories(
    kinova_utils/include
    kinova_system/include
    kinova_filter/include
    kinova_dynamics/include
    ${EIGEN3_INCLUDE_DIR}
    ${TINYXML2_INCLUDE_DIR}
)

# Sub directories
add_subdirectory(kinova_utils)
add_subdirectory(kinova_system)
add_subdirectory(kinova_filter)
add_subdirectory(kinova_dynamics)

# Add the executable
add_executable(roahm_kinova main.cpp)

# Link against libraries
target_link_libraries(roahm_kinova
    kinova_utils
    kinova_system
    kinova_filter
    kinova_dynamics
    Eigen3::Eigen
    fmt::fmt
    tinyxml2
)

# Set include directories for the target
target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_SOURCE_DIR}/kinova_utils/include
    ${CMAKE_SOURCE_DIR}/kinova_system/include
    ${CMAKE_SOURCE_DIR}/kinova_filter/include
    ${CMAKE_SOURCE_DIR}/kinova_dynamics/include
)