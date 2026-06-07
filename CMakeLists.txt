cmake_minimum_required(VERSION 3.19)

# stops CMake from trying to write files into the source directory
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed.")
endif()

set(NILT_VERSION "2.0.0")
project(nilt VERSION ${NILT_VERSION} LANGUAGES CXX)

# options
option(NILT_BUILD_EXAMPLES      "Build the example executables"       ON)
option(NILT_BUILD_TESTS         "Build and register CTest tests"      OFF)
option(NILT_BUILD_VERIFICATION  "Build the verification executables"  OFF)

# header-only library target
add_library(nilt INTERFACE)
add_library(nilt::nilt ALIAS nilt)

target_include_directories(nilt INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(nilt INTERFACE cxx_std_14)

add_library(nilt_warnings INTERFACE)
if(MSVC)
    target_compile_options(nilt_warnings INTERFACE /W4)
else()
    target_compile_options(nilt_warnings INTERFACE -Wall -Wextra -Wpedantic)
endif()

# tests
if(NILT_BUILD_TESTS)
    enable_testing()

    include(FetchContent)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.9.1
    )
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
    include(Catch)

    add_executable(nilt_tests
        tests/test_stehfest.cpp
        tests/test_talbot.cpp
        tests/test_dehoog.cpp
        tests/test_verification.cpp
    )
    target_link_libraries(nilt_tests PRIVATE nilt nilt_warnings Catch2::Catch2WithMain)
    target_include_directories(nilt_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    target_include_directories(nilt_tests SYSTEM PRIVATE ${catch2_SOURCE_DIR}/src)

    catch_discover_tests(nilt_tests)
endif()

# verification
if(NILT_BUILD_VERIFICATION)
    add_executable(verification          verification/verification.cpp)
    add_executable(benchmark_timing      verification/benchmark_timing.cpp)

    target_link_libraries(verification          PRIVATE nilt nilt_warnings)
    target_link_libraries(benchmark_timing      PRIVATE nilt nilt_warnings)

    target_include_directories(verification PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/verification"
    )
endif()

# examples
if(NILT_BUILD_EXAMPLES)
    add_executable(transport_sphere_diffusion    examples/transport/sphere_diffusion.cpp)
    add_executable(transport_cylinder_diffusion  examples/transport/cylinder_diffusion.cpp)
    add_executable(transport_advection_plume_2d  examples/transport/advection_plume_2d.cpp)
    add_executable(groundwater_theis_well        examples/groundwater/theis_well.cpp)
    add_executable(groundwater_well_dipole       examples/groundwater/well_dipole.cpp)

    target_link_libraries(transport_sphere_diffusion    PRIVATE nilt nilt_warnings)
    target_link_libraries(transport_cylinder_diffusion  PRIVATE nilt nilt_warnings)
    target_link_libraries(transport_advection_plume_2d  PRIVATE nilt nilt_warnings)
    target_link_libraries(groundwater_theis_well        PRIVATE nilt nilt_warnings)
    target_link_libraries(groundwater_well_dipole       PRIVATE nilt nilt_warnings)
endif()

# install rules
include(GNUInstallDirs)

install(TARGETS nilt EXPORT niltTargets)

install(DIRECTORY include/
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING PATTERN "*.hpp"
)

install(EXPORT niltTargets
    FILE        niltTargets.cmake
    NAMESPACE   nilt::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/nilt"
)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/niltConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/niltConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/niltConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/nilt"
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/niltConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/niltConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/nilt"
)
