include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(bench-fft fftw/bench-fft.cpp)
target_link_libraries(bench-fft argparse fftw-cpp)

add_executable(bench-fft-multi fftw/bench-fft-multi.cpp)
target_link_libraries(bench-fft-multi argparse fftw-cpp)

