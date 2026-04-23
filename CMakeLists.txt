cmake_minimum_required(VERSION 3.14)
set(NILT_VERSION "1.4.0")
project(nilt VERSION ${NILT_VERSION} LANGUAGES CXX)

# options
option(NILT_BUILD_EXAMPLES  "Build the example executables"  ON )
option(NILT_BUILD_TESTS     "Build and register CTest tests" OFF)

# header-only library target
add_library(nilt INTERFACE)
add_library(nilt::nilt ALIAS nilt)

target_include_directories(nilt INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(nilt INTERFACE cxx_std_14)

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
    target_link_libraries(nilt_tests PRIVATE nilt Catch2::Catch2WithMain)
    if(NOT MSVC)
        target_compile_options(nilt_tests PRIVATE -Wall -Wextra -Wpedantic)
    endif()

    catch_discover_tests(nilt_tests)
endif()

# examples
if(NILT_BUILD_EXAMPLES)
    #   target_name                   source_path
    set(_nilt_examples
        verification                  examples/verification/verification.cpp
        benchmark_timing              examples/verification/benchmark_timing.cpp
        transport_sphere_diffusion    examples/transport/sphere_diffusion.cpp
        transport_cylinder_diffusion  examples/transport/cylinder_diffusion.cpp
        transport_advection_plume_2d  examples/transport/advection_plume_2d.cpp
        groundwater_theis_well        examples/groundwater/theis_well.cpp
        groundwater_well_dipole       examples/groundwater/well_dipole.cpp
    )
    list(LENGTH _nilt_examples _len)
    math(EXPR _last "${_len} - 1")
    foreach(_i RANGE 0 ${_last} 2)
        math(EXPR _j "${_i} + 1")
        list(GET _nilt_examples ${_i} _target)
        list(GET _nilt_examples ${_j} _source)
        add_executable(${_target} ${_source})
        target_link_libraries(${_target} PRIVATE nilt)
        if(NOT MSVC)
            target_compile_options(${_target} PRIVATE -Wall -Wextra -Wpedantic)
        endif()

        # Place each binary next to its sources in a build/ subdirectory
        get_filename_component(_source_dir ${_source} DIRECTORY)
        set(_build_dir "${PROJECT_SOURCE_DIR}/${_source_dir}/build")
        set_target_properties(${_target} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${_build_dir}"
        )
        file(MAKE_DIRECTORY "${_build_dir}")
    endforeach()

    # verification needs the test headers in its own directory
    target_include_directories(verification PRIVATE
        "${PROJECT_SOURCE_DIR}/examples/verification"
    )
endif()

# install rules
include(GNUInstallDirs)

install(TARGETS nilt EXPORT niltTargets)

install(DIRECTORY include/
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
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
