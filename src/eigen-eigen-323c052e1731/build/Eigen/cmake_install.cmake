# Install script for directory: /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE FILE FILES
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/Cholesky"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/CholmodSupport"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/Core"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/Dense"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/Eigen"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/Eigenvalues"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/Geometry"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/Householder"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/IterativeLinearSolvers"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/Jacobi"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/LU"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/MetisSupport"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/OrderingMethods"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/PaStiXSupport"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/PardisoSupport"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/QR"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/QtAlignedMalloc"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/SPQRSupport"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/SVD"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/Sparse"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/SparseCholesky"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/SparseCore"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/SparseLU"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/SparseQR"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/StdDeque"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/StdList"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/StdVector"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/SuperLUSupport"
    "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/UmfPackSupport"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE DIRECTORY FILES "/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

