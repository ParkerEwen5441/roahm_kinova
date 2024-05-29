cmake_minimum_required(VERSION 3.10)

# Project name
project(kinova_filter)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(
    include
    include/kinova_filter
    ${CMAKE_SOURCE_DIR}/../kinova_system/include
)

# Add an interface library for the header-only library
add_library(${PROJECT_NAME} INTERFACE)

# Set include directories for the target
target_include_directories(${PROJECT_NAME} INTERFACE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/kinova_filter
    ${CMAKE_SOURCE_DIR}/../kinova_system/include
)