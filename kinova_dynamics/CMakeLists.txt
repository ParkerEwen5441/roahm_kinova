cmake_minimum_required(VERSION 3.10)

# Project name
project(kinova_dynamics)

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
    ${CMAKE_SOURCE_DIR}/kinova_dynamics/include
    ${CMAKE_SOURCE_DIR}/kinova_dynamics/include/kinova_dynamics
    ${CMAKE_SOURCE_DIR}/kinova_dynamics/include/kinova_dynamics/spatial
    ${CMAKE_SOURCE_DIR}/kinova_utils/include
    ${CMAKE_SOURCE_DIR}/kinova_system/include
    ${EIGEN3_INCLUDE_DIR}
    ${TINYXML2_INCLUDE_DIR}
)

# Source files
set(SOURCE_FILES
    src/KinovaControllerBlock.cpp
    src/ControllerBlock.cpp
    src/Helper.cpp
    src/Model.cpp
    src/Parser.cpp
    src/RNEABlock.cpp
    src/RobustControlBlock.cpp
    src/Spatial.cpp
)

# Add library
add_library(kinova_dynamics ${SOURCE_FILES})

# Install
install(TARGETS kinova_dynamics EXPORT KinovaDynamicsTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
)
export(EXPORT KinovaDynamicsTargets FILE KinovaDynamicsTargets.cmake)

# Link libraries
target_link_libraries(kinova_dynamics PUBLIC
    Eigen3::Eigen
    ${TINYXML_LIB}
    fmt::fmt
    ${CMAKE_SOURCE_DIR}/kinova_utils/build/libkinova_utils.a
    ${CMAKE_SOURCE_DIR}/kinova_system/build/libkinova_system.a
)

# Set include directories for the target
target_include_directories(kinova_dynamics PUBLIC
    ${CMAKE_SOURCE_DIR}/kinova_dynamics/include
    ${CMAKE_SOURCE_DIR}/kinova_dynamics/include/kinova_dynamics
    ${CMAKE_SOURCE_DIR}/kinova_dynamics/include/kinova_dynamics/spatial
    ${CMAKE_SOURCE_DIR}/kinova_utils/include
    ${CMAKE_SOURCE_DIR}/kinova_system/include
)