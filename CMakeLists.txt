cmake_minimum_required(VERSION 3.20)
project(ALGBIO_Lab_1)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(googletest)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

include_directories(BMSTU_NMLA_Lab_1)

add_executable(main.cpp main.cpp)

add_executable(tests.cpp tests.cpp)
target_link_libraries(tests.cpp gtest gtest_main)
