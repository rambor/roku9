cmake_minimum_required(VERSION 3.8)
set(roku9_LANGUAGES CXX)
set(roku9_VERSION 0.7)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -g1 -DNDEBUG -O3 -msse4.1")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -msse4.1")


#SET(BUILD_SHARED_LIBS ON)
project(roku9 LANGUAGES ${roku9_LANGUAGES} VERSION ${roku9_VERSION})
find_library(SASTools REQUIRED)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/src/cpd-master/build/cmake")
MESSAGE("MODULES PATH ${CMAKE_MODULE_PATH}")

set(Cpd_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/cpd-master/build/cmake")
find_package(Cpd REQUIRED)

SET (QHULL_INCLUDE_DIRS ~/src_libs/qhull-2015.2/src ~/src_libs/qhull-2015.2/src/libqhull ~/src_libs/qhull-2015.2/src/libqhullcpp)
include_directories(${QHULL_INCLUDE_DIRS})

find_library(QHULL_LIB libqhullstatic.a PATHS ~/src_libs/qhull-2015.2/lib)
find_library(QHULLCPP_LIB libqhullcpp.a PATHS ~/src_libs/qhull-2015.2/lib)


#
# BOOST SECTION
#
# set Boost directory
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_STATIC_RUNTIME ON)
# Boost root variable set as environment variable in CLION Preferences
MESSAGE("** BOOST ROOT => ${BOOSTROOT}")
SET(BOOST_ROOT ${BOOSTROOT})
SET(BOOST_LIBRARYDIR ${BOOSTROOT}/stage/lib)
SET(BOOST_INCLUDEDIR ${BOOSTROOT})
include_directories(${BOOST_INCLUDEDIR})

set(BOOST_FILESYSTEM ${BOOSTROOT}/stage/lib/libboost_filesystem.a )
set(BOOST_PO ${BOOSTROOT}/stage/lib/libboost_program_options.a )
set(BOOST_RANDOM ${BOOSTROOT}/stage/lib/libboost_random.a )
set(BOOST_REGEX ${BOOSTROOT}/stage/lib/libboost_regex.a )
set(BOOST_SYSTEM ${BOOSTROOT}/stage/lib/libboost_system.a )

MESSAGE("** BOOST INCLUDE => ${BOOST_INCLUDEDIR}")
MESSAGE("** BOOST SPECIAL => ${BOOST_FILESYSTEM}")
MESSAGE("** BOOST SPECIAL => ${BOOST_PO}")
MESSAGE("** BOOST SPECIAL => ${BOOST_RANDOM}")
MESSAGE("** BOOST SPECIAL => ${BOOST_REGEX}")
MESSAGE("** BOOST SPECIAL => ${BOOST_SYSTEM}")
MESSAGE("** BOOST INCLUDE DIRECTORY: ${Boost_INCLUDE_DIR}")

# when building, don't use the install RPATH already
# (but later on when installing)
MESSAGE("**  CMAKE INSTALL RPATH => ${CMAKE_INSTALL_RPATH}")
MESSAGE("** CMAKE INSTALL PREFIX => ${CMAKE_INSTALL_PREFIX}")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

MESSAGE("MODULES PATH ${CMAKE_MODULE_PATH}")

SET(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/eigen-eigen-323c052e1731)
include_directories(${EIGEN3_INCLUDE_DIR})

set(SOURCE_FILES main.cpp)
include_directories(src)
add_subdirectory(src)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_INSTALL_RPATH "../lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
MESSAGE("** CMAKE INSTALL RPATH => ${CMAKE_INSTALL_RPATH}")

#SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,../lib -l")
#MESSAGE("** CMAKE_EXE_LINKER_FLAGS => ${CMAKE_EXE_LINKER_FLAGS}")
set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "")

add_executable(roku9 ${SOURCE_FILES})
target_link_libraries(roku9
        PUBLIC
        Cpd::Library-C++
        ${BOOST_FILESYSTEM} ${BOOST_PO} ${BOOST_RANDOM} ${BOOST_REGEX} ${BOOST_SYSTEM}
        ${QHULL_LIB} ${QHULLCPP_LIB} Roku9 SASTools
        )

#set_target_properties(roku9 PROPERTIES LINK_OPTIONS "-Wl;-rpath;./;-nostdinc;-isystem;-cxx-isystem;-L../lib")
#set_target_properties(roku9 PROPERTIES INSTALL_RPATH "@loader_path/../lib")
#set_target_properties(roku9 PROPERTIES
#        ARCHIVE_OUTPUT_DIRECTORY
#        "${CMAKE_BINARY_DIR}/lib"
#        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
#        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
#        )

#link_directories(${CMAKE_BINARY_DIR}/lib)

MESSAGE("** CMAKE BIN DIR => ${CMAKE_BINARY_DIR}")

# on install, create distribution directory and bin and lib subdirectores
# copy dylib to lib and link binaries to it using rpath
enable_testing()
add_subdirectory("${PROJECT_SOURCE_DIR}/lib")

# this will install libs into directory
#install(TARGETS roku9 CONFIGURATIONS Release
#        RUNTIME DESTINATION "${CMAKE_BINARY_DIR}/bin"
#        LIBRARY DESTINATION "${CMAKE_BINARY_DIR}/lib"
#)

#add_custom_target(install_${PROJECT_NAME} "${CMAKE_COMMAND}" --build "${CMAKE_BINARY_DIR}" --target install
#        DEPENDS ${PROJECT_NAME}
#        COMMENT "Installing ${PROJECT_NAME}")