cmake_minimum_required(VERSION 3.10)

# Project name
project(kinova_system)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(
    include
    include/kinova_system
    include/kinova_system/Port
    ${CMAKE_SOURCE_DIR}/../kinova_utils/include
)

# Source files
set(SOURCE_FILES
    src/BaseBlock.cpp
    src/System.cpp
)

# Add library
add_library(${PROJECT_NAME} ${SOURCE_FILES})

# Link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC rclcpp)

# Set include directories for the target
target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/kinova_system
    ${CMAKE_SOURCE_DIR}/include/kinova_system/Port
    ${CMAKE_SOURCE_DIR}/../kinova_utils/include
)