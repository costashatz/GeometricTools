cmake_minimum_required (VERSION 2.6)
project (LinearAlgebraTools)

set(ROOT .)

set(CMAKE_CXX_FLAGS "-std=c++11")

include_directories(${ROOT}/includes)

add_subdirectory(Example)