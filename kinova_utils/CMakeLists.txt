cmake_minimum_required(VERSION 3.10)

# Project name and version
project(kinova_utils VERSION 1.0 LANGUAGES CXX)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include directories
include_directories(include)

# Find the required packages
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(spdlog REQUIRED)

# Source files
set(SOURCES
    src/Utils.cpp
)

# Header files
set(HEADERS
    include/kinova_utils/Logger.hpp
    include/kinova_utils/Utils.hpp
)

# Create the library
add_library(kinova_utils STATIC ${SOURCES} ${HEADERS})

# Link the library with the required packages
target_link_libraries(kinova_utils PUBLIC Eigen3::Eigen spdlog::spdlog)