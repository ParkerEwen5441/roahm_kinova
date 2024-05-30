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
find_package(fmt REQUIRED)

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/kinova_utils/include
    ${CMAKE_SOURCE_DIR}/kinova_utils/include/kinova_utils
)

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
target_include_directories(kinova_utils PUBLIC ${CMAKE_SOURCE_DIR}/kinova_utils/include)

# Install
install(TARGETS kinova_utils EXPORT KinovaUtilsTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
)
export(EXPORT KinovaUtilsTargets FILE KinovaUtilsTargets.cmake)


# Link the library with the required packages
target_link_libraries(kinova_utils PUBLIC Eigen3::Eigen spdlog::spdlog fmt::fmt)

# Set include directories for the target
target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_SOURCE_DIR}/kinova_utils/include
    ${CMAKE_SOURCE_DIR}/kinova_utils/include/kinova_utils
)