cmake_minimum_required(VERSION 3.10)
project(nilt VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(benchmark module/benchmark.cpp)
target_include_directories(benchmark PUBLIC "${PROJECT_SOURCE_DIR}"        )
target_include_directories(benchmark PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_include_directories(benchmark PUBLIC "${PROJECT_SOURCE_DIR}/tests"  )
