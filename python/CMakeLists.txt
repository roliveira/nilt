cmake_minimum_required(VERSION 3.14)
project(nilt_python LANGUAGES CXX)

find_package(pybind11 REQUIRED)

pybind11_add_module(_nilt nilt_bind.cpp)

target_include_directories(_nilt PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../include"
)

target_compile_features(_nilt PRIVATE cxx_std_14)

install(TARGETS _nilt DESTINATION nilt)
