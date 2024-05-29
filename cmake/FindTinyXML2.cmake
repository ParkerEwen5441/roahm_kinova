# Locate the TinyXML2 library

find_path(TINYXML2_INCLUDE_DIR NAMES tinyxml2.h PATH_SUFFIXES tinyxml2)
find_library(TINYXML2_LIBRARY NAMES tinyxml2)

# Report the result
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TinyXML2 REQUIRED_VARS TINYXML2_INCLUDE_DIR TINYXML2_LIBRARY)

