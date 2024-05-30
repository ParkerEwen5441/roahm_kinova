cmake_minimum_required(VERSION 3.10)

# Project name
project(kinova_system)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(fmt REQUIRED)

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/kinova_system/include
    ${CMAKE_SOURCE_DIR}/kinova_system/include/kinova_system
    ${CMAKE_SOURCE_DIR}/kinova_system/include/kinova_system/Port
    ${CMAKE_SOURCE_DIR}/kinova_utils/include
)

# Source files
set(SOURCE_FILES
    src/BaseBlock.cpp
    src/System.cpp
)

# Add library
add_library(${PROJECT_NAME} ${SOURCE_FILES})

# Install
install(TARGETS kinova_system EXPORT KinovaSystemTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
)
export(EXPORT KinovaSystemTargets FILE KinovaSystemTargets.cmake)

# Link against libraries
target_link_libraries(kinova_system
    fmt::fmt
)

# Set include directories for the target
target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_SOURCE_DIR}/kinova_system/include
    ${CMAKE_SOURCE_DIR}/kinova_system/include/kinova_system
    ${CMAKE_SOURCE_DIR}/kinova_system/include/kinova_system/Port
    ${CMAKE_SOURCE_DIR}/kinova_utils/include
)