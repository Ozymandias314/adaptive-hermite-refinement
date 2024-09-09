cmake_minimum_required(VERSION 3.21)
set(CMAKE_SYSTEM_PROCESSOR arm64-apple-darwin) # TODO: not for everything
project(adaptive_hermite_refinement)

set(CMAKE_CXX_STANDARD 20)

# ==========================
# Library dependencies
# TODO: find locally first if available
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
        GIT_TAG 9ceface91483775a6c74d06ebf717bbb2768452f
)
FetchContent_MakeAvailable(mdspan)

FetchContent_Declare(
        cnpy
        GIT_REPOSITORY https://github.com/rogersce/cnpy.git
        GIT_TAG 4e8810b1a8637695171ed346ce68f6984e585ef4
)
FetchContent_MakeAvailable(cnpy)

target_include_directories(cnpy PUBLIC ${CNPY_SOURCE_DIR})

# ==========================
# Project code
# ==========================
enable_testing()
add_subdirectory(fftw-cpp)
target_link_libraries(fftw-cpp INTERFACE mdspan)

option(ENABLE_AVX512 "Enable AVX512 instructions (e.g. might want to turn off for the sake of Valgrind)" ON)
if (ENABLE_AVX512)
    add_compile_options("-mavx512f" "-march=native" "-ffast-math")
endif ()
add_compile_options("-fopencilk")
add_link_options("-fopencilk")

option(CILK_STUB "Use stubs instead of Cilk constructs" OFF)
if (CILK_STUB)
    add_compile_definitions(CILK_STUB)
endif ()

add_executable(adaptive_hermite_refinement test-triangle.cpp)
target_link_libraries(adaptive_hermite_refinement argparse mdspan)

add_executable(bench-fft bench-fft.cpp)
target_link_libraries(bench-fft argparse fftw-cpp)

add_executable(bench-fft-multi bench-fft-multi.cpp)
target_link_libraries(bench-fft-multi argparse fftw-cpp)

add_library(src-lib OBJECT Naive.cpp Naive.h HermiteRunner.cpp HermiteRunner.h)
target_include_directories(src-lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(src-lib mdspan fftw-cpp cnpy)

add_executable(naive driver.cpp)
target_link_libraries(naive argparse src-lib)

option(AHR_TEST "Enable AHR tests" ON)
if (AHR_TEST)
    enable_testing()
    add_subdirectory(test)
endif ()
