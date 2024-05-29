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

# Include directories
include_directories(
    include
    include/kinova_dynamics
    include/kinova_dynamics/spatial
    ${CMAKE_SOURCE_DIR}/../kinova_utils/include
    ${CMAKE_SOURCE_DIR}/../kinova_system/include
    ${EIGEN3_INCLUDE_DIR}
    ${TINYXML2_INCLUDE_DIR}
)

# Source files
set(SOURCE_FILES
    src/ControllerBlock.cpp
    src/Helper.cpp
    src/Model.cpp
    src/Parser.cpp
    src/RNEABlock.cpp
    src/RobustControlBlock.cpp
    src/Spatial.cpp
)

# Add library
add_library(${PROJECT_NAME} ${SOURCE_FILES})

# Link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC
    Eigen3::Eigen
    ${TINYXML_LIB}
    ${CMAKE_SOURCE_DIR}/../kinova_utils/libkinova_utils.a
    ${CMAKE_SOURCE_DIR}/../kinova_system/libkinova_system.a
)

# Set include directories for the target
target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/kinova_dynamics
    ${CMAKE_SOURCE_DIR}/include/kinova_dynamics/spatial
    ${CMAKE_SOURCE_DIR}/../kinova_utils/include
    ${CMAKE_SOURCE_DIR}/../kinova_system/include
)