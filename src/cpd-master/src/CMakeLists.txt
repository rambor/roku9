cmake_minimum_required(VERSION 3.8)
project(src)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        affine.cpp
        gauss_transform.cpp
        gauss_transform_fgt.cpp
        gauss_transform_make_default.cpp
        matrix.cpp
        nonrigid.cpp
        normalization.cpp
        rigid.cpp
        transform.cpp
        utils.cpp)

add_executable(src ${SOURCE_FILES})