cmake_minimum_required(VERSION 3.14)
project(nilt_python LANGUAGES CXX)

# Read the canonical version from the parent CMakeLists.txt so the
# pybind11 module and the C++ library never disagree.
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/../CMakeLists.txt" _nilt_root_cmake)
string(REGEX MATCH "NILT_VERSION[ \t]+\"([0-9]+)\\.([0-9]+)\\.([0-9]+)\""
       _nilt_ver_match "${_nilt_root_cmake}")
if(NOT _nilt_ver_match)
    message(FATAL_ERROR "Could not extract NILT_VERSION from ../CMakeLists.txt")
endif()
set(NILT_VERSION_MAJOR "${CMAKE_MATCH_1}")
set(NILT_VERSION_MINOR "${CMAKE_MATCH_2}")
set(NILT_VERSION_PATCH "${CMAKE_MATCH_3}")
set(NILT_VERSION "${NILT_VERSION_MAJOR}.${NILT_VERSION_MINOR}.${NILT_VERSION_PATCH}")

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/../include/nilt/version.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/generated/nilt/version.hpp"
    @ONLY
)

find_package(pybind11 REQUIRED)

pybind11_add_module(_nilt nilt_bind.cpp)

target_include_directories(_nilt PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../include"
    "${CMAKE_CURRENT_BINARY_DIR}/generated"
)

target_compile_features(_nilt PRIVATE cxx_std_14)

install(TARGETS _nilt DESTINATION nilt)
