include(FetchContent)

FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG 1d17ea141d2c11b8917d2c7d029f1c4e2b9769b2
)
FetchContent_MakeAvailable(googletest)

add_executable(unit-tests naive.cpp)
target_link_libraries(unit-tests src-lib GTest::gmock_main)
