#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Cpd::Library-C++" for configuration ""
set_property(TARGET Cpd::Library-C++ APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(Cpd::Library-C++ PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libcpd.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS Cpd::Library-C++ )
list(APPEND _IMPORT_CHECK_FILES_FOR_Cpd::Library-C++ "${CMAKE_CURRENT_SOURCE_DIR}/src/cpd-master/build/libcpd.a" )
# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
