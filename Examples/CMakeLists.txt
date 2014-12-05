cmake_minimum_required (VERSION 2.6)
project (LinearAlgebraTools)


add_executable(Example main.cpp)
target_link_libraries(Example ${PROJECT_NAME})
