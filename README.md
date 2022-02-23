## Geometric Tools

### Math Library for Linear Algebra, Numerical Methods, Geometric Representations and Collision Detection Algorithms

#### Dependencies
1. **C++11 compiler**
	* Tested on gcc>=4.6, clang>=3.1
	* Should compile on Visual Studio 2013-2015
2. **CMake**
	* Tested on CMake>=2.8
	* Should work on 2.6 too
3. **LAPACK/BLAS**
	* Needed for Linear Algebra operations
	* In *ubuntu Linux, just install the following packages: liblapack-dev libblas-dev
4. **Google Tests (gtest)** *[OPTIONAL]*
	* Needed to build the test codes

#### How to compile:

```
cd path/to/GeometricTools
mkdir build && cd build
cmake ..
make
```

* Optional CMake options:
	1. BUILD_TEST (ON/OFF) - Determine to either build tests or not. Defaults to OFF.
	2. BUILD_EXAMPLES (ON/OFF) - Specify whether you want to build tests or not. Defaults to OFF.
	3. RUN_TEST (ON/OFF) - Specify whether you want to automatically execute all tests upon build. Defaults to OFF. 

#### How to use:

* You can just edit the CMakeLists.txt and include the headers to your project (after all it is just a header library).
* Via ExternalProject to always get latest edition:
	```cmake
	find_package(LAPACK REQUIRED)
	set( ENV{BLA_VENDOR} "Generic" )
	find_package(BLAS REQUIRED)

	include(ExternalProject)
	ExternalProject_Add(GeometricToolsProj
	    GIT_REPOSITORY "https://github.com/costashatz/GeometricTools"
	    CMAKE_ARGS -DBUILD_TEST=OFF -DBUILD_EXAMPLES=OFF
	    PREFIX "${CMAKE_CURRENT_BINARY_DIR}"
	    INSTALL_COMMAND ""
	)
	# Specify include dir
	ExternalProject_Get_Property(GeometricToolsProj source_dir)
	set(GEOMETRIC_TOOLS_INCLUDE_DIRS ${source_dir}/include)

	# Specify MainTest's link libraries
	ExternalProject_Get_Property(GeometricToolsProj binary_dir)
	set(GEOMETRIC_TOOLS_LIBS_DIR ${binary_dir})

	link_directories(${GEOMETRIC_TOOLS_LIBS_DIR})
	include_directories(${PROJECT_SOURCE_DIR}/include ${GEOMETRIC_TOOLS_INCLUDE_DIRS})
	....
	target_link_libraries(mytarget GeometricTools ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
	```


#### So far I have implemented:

1. Generic Templated Class for Matrices
    * Supports any dimension and many properties
2. Generic Templated Class for Vectors
    * Supports any dimension and most of the vector properties
3. Solve Linear Systems
    * Using Gauss Elimination - **error prone**
    * Using LU Decomposition
4. Decompositions
    * LU Decomposition
    * QR Decomposition
    * SVD Decomposition - **not working right now**
5. Linear Shapes
	* Classes for basic linear shapes (line, ray, segment)
6. Polygons
	* Classes for basic 2D polygons (triangle, rectangle, polyline, general polygons)
7. Curves
	* Specific quadratic curves (defined by xTAx+bTx+c=0)
	* Generic polynomial curves/splines (templated in size [biggest power of curve]) - 1D functions
	* Piecewise Curves
	* Hermite Cubic Curves
	* Cardinal Cubic Curves
	* Plane (2D) Curves
8. Distances
	* Functions for computing distances between primitives (point-linear shapes, point-polyline, linear shapes-linear shapes, segment-polyline)
9. Intersections
	* Functions for computing intersection points (2D) between primitives (linear shapes-linear shapes, linear shapes-polygon)
10. Simple gtests (more need to be added)
11. Simple examples (more interactive ones need to be added)

#### Planning to implement:

1. More primitives (both in 2D and 3D)
	* Piecewise polynomial curves
	* Space(3D) polynomial curves
2. More Distances, Intersections in 2D/3D
3. Cleanify code/add comments


