# - Find GMP
# The GNU Multiple Precision Arithmetic Library
#
# The module defines the following variables:
#  GMP_FOUND - the system has GMP
#  GMP_INCLUDE_DIR - where to find gmp.h
#  GMP_INCLUDE_DIRS - gmp includes
#  GMP_LIBRARY - where to find the GMP library
#  GMP_LIBRARIES - aditional libraries
#  GMP_ROOT_DIR - root dir (ex. /usr/local)

find_path (GMP_INCLUDE_DIR NAMES gmp.h DOC "GMP include directory")

set (GMP_INCLUDE_DIRS ${GMP_INCLUDE_DIR})

find_library (GMP_LIBRARY NAMES gmp DOC "GMP library location")

set (GMP_LIBRARIES ${GMP_LIBRARY})

# root dir
# try to guess root dir from include dir
if (GMP_INCLUDE_DIR)
  string (REGEX REPLACE "(.*)/include.*" "\\1" GMP_ROOT_DIR ${GMP_INCLUDE_DIR})

# try to guess root dir from library dir
elseif (GMP_LIBRARY)
  string (REGEX REPLACE "(.*)/lib[/|32|64].*" "\\1" GMP_ROOT_DIR ${GMP_LIBRARY})
endif ()

# handle REQUIRED and QUIET options
include (FindPackageHandleStandardArgs)


find_package_handle_standard_args (GMP
  REQUIRED_VARS GMP_LIBRARY GMP_INCLUDE_DIR GMP_INCLUDE_DIRS GMP_LIBRARIES
)


mark_as_advanced (
  GMP_LIBRARY
  GMP_LIBRARIES
  GMP_INCLUDE_DIR
  GMP_INCLUDE_DIRS
  GMP_ROOT_DIR
)
