# Install script for directory: /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/AdolcForward"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/AlignedVector3"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/ArpackSupport"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/AutoDiff"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/BVH"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/EulerAngles"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/FFT"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/IterativeSolvers"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/KroneckerProduct"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/LevenbergMarquardt"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/MatrixFunctions"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/MoreVectorization"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/MPRealSupport"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/NonLinearOptimization"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/NumericalDiff"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/OpenGLSupport"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/Polynomials"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/Skyline"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/SparseExtra"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/SpecialFunctions"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/Splines"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

