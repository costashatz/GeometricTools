cmake_minimum_required (VERSION 2.6)
project (GeometricTools)


add_executable(LinearSystems main.cpp)
target_link_libraries(LinearSystems ${PROJECT_NAME})
