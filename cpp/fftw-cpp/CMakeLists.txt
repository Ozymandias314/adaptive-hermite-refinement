project(fftw-cpp LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)

add_library(fftw-cpp INTERFACE)
target_link_libraries(fftw-cpp INTERFACE fftw3 m)
target_include_directories(fftw-cpp INTERFACE include/)

function(fftw_test NAME)
add_executable(${NAME} ${ARGN})
target_link_libraries(${NAME} fftw-cpp)
add_test(${NAME} ${NAME})
endfunction()

fftw_test(test-fftw test.cpp)
fftw_test(test-fftw-wrapper test-wrapper.cpp)
fftw_test(test-fftw-wrapper-2d test-wrapper-2d.cpp)
fftw_test(test-fftw-r2c-2d test-r2c-2d.cpp)
fftw_test(test-plan test-plan.cpp)
fftw_test(test-mdspan test-mdspan.cpp)
