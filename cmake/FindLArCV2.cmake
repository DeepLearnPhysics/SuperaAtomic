# The baseline LArCV2 UPS package does not export a CMake target.
find_path(LArCV2_INCLUDE_DIR larcv/core/DataFormat/Flash.h HINTS ENV LARCV2_INC)
find_library(LArCV2_LIBRARY NAMES larcv HINTS ENV LARCV2_LIB)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LArCV2 REQUIRED_VARS LArCV2_INCLUDE_DIR LArCV2_LIBRARY)
if(NOT TARGET LArCV2::larcv)
  add_library(LArCV2::larcv SHARED IMPORTED)
  set_target_properties(LArCV2::larcv PROPERTIES
    IMPORTED_LOCATION "${LArCV2_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${LArCV2_INCLUDE_DIR}")
endif()
mark_as_advanced(LArCV2_INCLUDE_DIR LArCV2_LIBRARY)
