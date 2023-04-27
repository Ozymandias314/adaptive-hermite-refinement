project(fftw-cpp LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

add_library(fftw-cpp INTERFACE)
target_link_libraries(fftw-cpp INTERFACE fftw3 m)
target_include_directories(fftw-cpp INTERFACE include/)

add_executable(test-fftw test.cpp)
target_link_libraries(test-fftw fftw-cpp)

add_executable(test-fftw-wrapper test-wrapper.cpp)
target_link_libraries(test-fftw-wrapper fftw-cpp)

add_executable(test-fftw-wrapper-2d test-wrapper-2d.cpp)
target_link_libraries(test-fftw-wrapper-2d fftw-cpp)
