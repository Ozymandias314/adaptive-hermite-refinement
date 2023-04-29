cmake_minimum_required(VERSION 3.21)
project(adaptive_hermite_refinement)

set(CMAKE_CXX_STANDARD 20)

# ==========================
# Library dependencies
# ==========================
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

# ==========================
# Project code
# ==========================
add_subdirectory(fftw-cpp)
target_link_libraries(fftw-cpp INTERFACE mdspan)

add_compile_options("-mavx512f" "-march=native" "-fopencilk" "-ffast-math")
add_link_options("-fopencilk")

option(CILK_STUB "Use stubs instead of Cilk constructs" OFF)
if(CILK_STUB)
    add_compile_definitions(CILK_STUB)
endif()

add_executable(adaptive_hermite_refinement test-triangle.cpp)
target_link_libraries(adaptive_hermite_refinement argparse mdspan)

add_executable(bench-fft bench-fft.cpp)
target_link_libraries(bench-fft argparse fftw-cpp)

add_executable(bench-fft-multi bench-fft-multi.cpp)
target_link_libraries(bench-fft-multi argparse fftw-cpp)

add_executable(naive driver.cpp Naive.cpp Naive.h HermiteRunner.cpp HermiteRunner.h)
target_link_libraries(naive argparse fftw-cpp)
