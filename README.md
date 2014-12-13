## Geometric Tools

###Math Library for Linear Algebra, Numerical Methods, Geometric Representations and Collision Detection Algorithms

####Needs C++11 compatible compiler to compile!

####How to compile:

```
cd path/to/LinearAlgebraTools
mkdir build && cd build
cmake ..
make
```

* Optional CMake options:
	1. BUILD_TEST (ON/OFF) - Determine to either build tests or not. Defaults to OFF.
	2. BUILD_EXAMPLES (ON/OFF) - Specify whether you want to build tests or not. Defaults to OFF.

####How to use:

* You can just edit the CMakeLists.txt and include the headers to your project (after all it is just a header library).
* Via ExternalProject to always get latest edition:
	```cmake
	include(ExternalProject)
	ExternalProject_Add(LinearAlgebraToolsProj
	    GIT_REPOSITORY "https://github.com/costashatz/LinearAlgebraTools"
	    CMAKE_ARGS -DBUILD_TEST=OFF -DBUILD_EXAMPLES=OFF
	    PREFIX "${CMAKE_CURRENT_BINARY_DIR}"
	    INSTALL_COMMAND ""
	)
	# Specify include dir
	ExternalProject_Get_Property(LinearAlgebraToolsProj source_dir)
	set(LINEAR_ALGEBRA_TOOLS_INCLUDE_DIRS ${source_dir}/include)

	# Specify MainTest's link libraries
	ExternalProject_Get_Property(LinearAlgebraToolsProj binary_dir)
	set(LINEAR_ALGEBRA_TOOLS_LIBS_DIR ${binary_dir})

	link_directories(${LINEAR_ALGEBRA_TOOLS_LIBS_DIR})
	include_directories(${PROJECT_SOURCE_DIR}/include ${LINEAR_ALGEBRA_TOOLS_INCLUDE_DIRS})
	....
	target_link_libraries(mytarget LinearAlgebraTools)
	```


####So far I have implemented:

1. Generic Templated Class for Matrices
    * Supports any dimension and many properties (many need to be added and some to improve)
2. Generic Templated Class for Vectors
    * Supports any dimension and most of the vector properties
3. Solve Linear Systems
    * Using Gauss Elimination
    * Using LU Decomposition
4. Decompositions
    * LU Decomposition
    * QR Decomposition
    * SVD Decomposition
5. Linear Shapes
	* Classes for basic 2D linear shapes (line, ray, segment)
6. Polygons
	* Classes for basic 2D polygons (triangle, rectangle, polyline, general polygons)
7. Curves
	* Specific quadratic curves (defined by xTAx+bTx+c=0)
	* Generic polynomial curves/splines (templated in size [biggest power of curve])
	* Piecewise Curves
	* Hermite Cubic Curves
	* Cardinal Cubic Curves
8. Distances in 2D
	* Functions for computing distances between 2D primitives (point-linear shapes, point-polyline, linear shapes-linear shapes, segment-polyline)

####Planning to implement:

1. More primitives (both in 2D and 3D)
	* Piecewise polynomial curves
	* Space(3D) polynomial curves
2. Distances, Intersections in 2D/3D
3. gtests
4. Improve the templated classes for Matrices and Vectors
    * Cleanify code/add comments
5. Add more linear systems solving mechanisms (algorithms)


####WORK IN PROGRESS - WILL BE UPDATED FREQUENTLY

I am currently working on this little library/project and will be pushing more updated stuff. Any suggestions of course are more than welcomed.