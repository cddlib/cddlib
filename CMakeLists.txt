cmake_minimum_required (VERSION 3.18)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

project (cddlib C)
list (APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})

option (BUILD_SHARED_LIBS "shared/static" ON)

include (GNUInstallDirs)

set(cdd_sources
  lib-src/cddcore.c
  lib-src/cddio.c
  lib-src/cddlib.c
  lib-src/cddlp.c
  lib-src/cddmp.c
  lib-src/cddproj.c
  lib-src/setoper.c)

add_library(cdd ${cdd_sources})

target_include_directories(cdd PRIVATE lib-src)


set_property(TARGET cdd PROPERTY POSITION_INDEPENDENT_CODE ON)
set_property(TARGET cdd PROPERTY WINDOWS_EXPORT_ALL_SYMBOLS ON)

install (TARGETS cdd
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install (FILES
  lib-src/cdd.h
  lib-src/cddmp.h
  lib-src/cddtypes.h
  lib-src/setoper.h
  lib-src/splitmix64.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cddlib)

  
function(generate_gmp_variant INPUT_FILE)
    get_filename_component(_name "${INPUT_FILE}" NAME)
    get_filename_component(_name_we "${INPUT_FILE}" NAME_WE)
    get_filename_component(_ext "${INPUT_FILE}" EXT)
    set(OUTPUT_FILE "${CMAKE_CURRENT_BINARY_DIR}/${_name_we}_f${_ext}")
    file(READ "${INPUT_FILE}" CONTENTS)
    string(REGEX REPLACE "dd_" "ddf_" CONTENTS "${CONTENTS}")
    string(REGEX REPLACE "cddf_" "cdd_" CONTENTS "${CONTENTS}")
    string(REGEX REPLACE "mytype" "myfloat" CONTENTS "${CONTENTS}")
    string(REGEX REPLACE "#include \"cdd.h\"" "#include \"cdd_f.h\"" CONTENTS "${CONTENTS}")
    string(REGEX REPLACE "#include \"cddtypes.h\"" "#include \"cddtypes_f.h\"" CONTENTS "${CONTENTS}")
    string(REGEX REPLACE "#include \"cddmp.h\"" "#include \"cddmp_f.h\"" CONTENTS "${CONTENTS}")
    string(REGEX REPLACE "__CDD_H" "__CDD_HF" CONTENTS "${CONTENTS}")
    string(REGEX REPLACE "__CDD_HFF" "__CDD_HF" CONTENTS "${CONTENTS}")
    string(REGEX REPLACE "__CDDMP_H" "_CDDMP_HF" CONTENTS "${CONTENTS}")
    string(REGEX REPLACE "__CDDTYPES_H" "_CDDTYPES_HF" CONTENTS "${CONTENTS}")
    string(REGEX REPLACE "GMPRATIONAL" "ddf_GMPRATIONAL" CONTENTS "${CONTENTS}")
    string(REGEX REPLACE "ARITHMETIC" "ddf_ARITHMETIC" CONTENTS "${CONTENTS}")
    string(REGEX REPLACE "CDOUBLE" "ddf_CDOUBLE" CONTENTS "${CONTENTS}")
    file(WRITE "${OUTPUT_FILE}" "/* generated automatically from ${INPUT_FILE} */\n${CONTENTS}")
    set(${ARGV1}_OUTPUT "${OUTPUT_FILE}" PARENT_SCOPE)
endfunction()


option (ENABLE_GMP "Enable GMP" ON)
if (ENABLE_GMP)
  find_package(GMP)
endif ()
if (GMP_FOUND)
  generate_gmp_variant (lib-src/cddcore.c)
  generate_gmp_variant (lib-src/cddlp.c)
  generate_gmp_variant (lib-src/cddmp.c)
  generate_gmp_variant (lib-src/cddio.c)
  generate_gmp_variant (lib-src/cddlib.c)
  generate_gmp_variant (lib-src/cddproj.c)
  generate_gmp_variant (lib-src/cdd.h)
  generate_gmp_variant (lib-src/cddmp.h)
  generate_gmp_variant (lib-src/cddtypes.h)

  add_library(cddgmp
    ${CMAKE_CURRENT_BINARY_DIR}/cddcore_f.c
    ${CMAKE_CURRENT_BINARY_DIR}/cddio_f.c
    ${CMAKE_CURRENT_BINARY_DIR}/cddlib_f.c
    ${CMAKE_CURRENT_BINARY_DIR}/cddlp_f.c
    ${CMAKE_CURRENT_BINARY_DIR}/cddmp_f.c
    ${CMAKE_CURRENT_BINARY_DIR}/cddproj_f.c
    ${cdd_sources})
  target_include_directories(cddgmp PRIVATE lib-src ${CMAKE_CURRENT_BINARY_DIR})
  target_compile_definitions(cddgmp PUBLIC GMPRATIONAL)
  set_property(TARGET cddgmp PROPERTY POSITION_INDEPENDENT_CODE ON)
  set_property(TARGET cddgmp PROPERTY WINDOWS_EXPORT_ALL_SYMBOLS ON)
  target_link_libraries(cddgmp PRIVATE ${GMP_LIBRARIES})

  install (TARGETS cddgmp
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  install (FILES
    ${CMAKE_CURRENT_BINARY_DIR}/cdd_f.h
    ${CMAKE_CURRENT_BINARY_DIR}/cddmp_f.h
    ${CMAKE_CURRENT_BINARY_DIR}/cddtypes_f.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cddlib)
endif ()

# pkgconfig
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix "\${prefix}")
set(libdir "\${prefix}/${CMAKE_INSTALL_LIBDIR}")
set(includedir "\${prefix}/${CMAKE_INSTALL_INCLUDEDIR}")
set(PACKAGE_NAME "cddlib")
set(PACKAGE_VERSION "0.94n")
configure_file(cddlib.pc.in cddlib.pc)
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/cddlib.pc
	 DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

# testing
add_executable(testcdd2 src/testcdd2.c)
target_include_directories(testcdd2 PRIVATE lib-src)
target_link_libraries(testcdd2 PRIVATE cdd)
enable_testing()
add_test(NAME testcdd2 COMMAND testcdd2)
