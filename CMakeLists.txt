cmake_minimum_required (VERSION 2.8.12)
project (GeometricTools)

include(CheckCXXCompilerFlag)

CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
	message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")

option(BUILD_TEST "Use Gtest to create the test cases for the code" OFF)
option(BUILD_EXAMPLES "Build examples of the code" OFF)

if(BUILD_TEST)
    # Setup testing
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/gtest)
    enable_testing()
    set(PROJECT_TEST_NAME ${PROJECT_NAME}_test)
    include_directories(${GTEST_INCLUDE_DIRS})
endif()

file(GLOB_RECURSE _srcs
    include/*.h
    include/*.hpp)

find_package( LAPACK REQUIRED)
include_directories(${LAPACK_INCLUDE_DIR})

add_library(${PROJECT_NAME} SHARED ${_srcs})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME} ${LAPACK_LIBRARIES})
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

if(BUILD_TEST)
    ##############
    # Unit Tests #
    ##############
    file(GLOB_RECURSE _test_srcs ${PROJECT_SOURCE_DIR}/test/*.cpp)
    add_executable(${PROJECT_TEST_NAME} ${_test_srcs})
    add_dependencies(${PROJECT_TEST_NAME} googletest)

    #target_link_libraries(${PROJECT_TEST_NAME} ${PROJECT_NAME})
    target_link_libraries(${PROJECT_TEST_NAME}
                          ${PROJECT_NAME}
                          ${GTEST_LIBS_DIR}/libgtest.a
                          ${GTEST_LIBS_DIR}/libgtest_main.a
                          pthread)

    # This is so you can do 'make test' to see all your tests run, instead of
    # manually running the executable runUnitTests to see those specific tests.
    add_test(NAME test
             COMMAND ${PROJECT_TEST_NAME})
endif()

if(BUILD_EXAMPLES)
  add_subdirectory(Examples)
endif()
