cmake_minimum_required(VERSION 3.22)
project(adaptive_hermite_refinement)

set(CMAKE_CXX_STANDARD 20)

include(FetchContent)
FetchContent_Declare(
        argparse
        GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

FetchContent_Declare(
        mdspan
        GIT_REPOSITORY https://github.com/kokkos/mdspan.git
        GIT_TAG        9ceface91483775a6c74d06ebf717bbb2768452f
)
FetchContent_MakeAvailable(mdspan)

add_executable(adaptive_hermite_refinement
        test-triangle.cpp)

target_link_libraries(adaptive_hermite_refinement argparse mdspan)
