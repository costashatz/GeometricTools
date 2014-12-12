cmake_minimum_required (VERSION 2.6)
project (LinearAlgebraTools)


add_executable(PlaneCurves main.cpp)
target_link_libraries(PlaneCurves ${PROJECT_NAME})
