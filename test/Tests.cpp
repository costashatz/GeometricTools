/**
Copyright (c) 2014, Konstantinos Chatzilygeroudis
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
    in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**/

#include <gtest/gtest.h>
#include <Math/Vector.h>
#include <Math/Matrix.h>
#include <Math/Transformations/2D/Homogeneous.h>
#include <Math/Transformations/3D/Homogeneous.h>
#include <Math/Numerical Optimization/1D/GoldenSearchMinimization.h>
#include <Primitives/LinearShapes.h>
#include <Primitives/Polyline.h>
#include <Primitives/2D/Triangle.h>
#include <Primitives/2D/Rectangle.h>
#include <Primitives/2D/Circle.h>
#include <Distances/PointToLinear.h>
#include <Distances/LinearToLinear.h>
#include <Distances/LinearToPolyline.h>

TEST(MathTest, VectorTests)
{
    using GeometricTools::Math::Vector;
    Vector<2> a{2,3}, b{5,6};
    EXPECT_EQ(a*b, 28.0);
    EXPECT_EQ(a+b, Vector<2>(7,9));
    EXPECT_EQ(a-b, Vector<2>(-3,-3));
    Vector<3> c{1,2,3}, d{4,5,6};
    EXPECT_EQ(c*d, 32.0);
    EXPECT_EQ(1.0/c, Vector<3>(1,0.5,1/3.0));
    EXPECT_DOUBLE_EQ(c.normalized().length(), 1.0);
    d.normalize();
    EXPECT_DOUBLE_EQ(d.length(), 1.0);
}

TEST(MathTet, MatrixTests)
{
    using GeometricTools::Math::Vector;
    using GeometricTools::Math::Matrix;
    typedef Matrix<2,2> Matrix2x2;
    Matrix2x2 A{1,0,0,1}, B{0,1,0,0}, ID = {1,0,0,1};
    Vector<2> c{4,5};
    EXPECT_EQ(A*c, Vector<2>(4,5));
    A = {1,1,0,1};
    EXPECT_EQ(A*c, Vector<2>(9,5));
    EXPECT_EQ(A-B, ID);
    EXPECT_EQ(A+B, Matrix2x2(1,2,0,1));
    EXPECT_EQ(A*B, Matrix2x2(0,1,0,0));
    EXPECT_DOUBLE_EQ(A.norm(), sqrt(3));
}

TEST(MathTest, Transformations2D)
{
    using GeometricTools::Math::Vector;
    using GeometricTools::Math::Matrix;
    using namespace GeometricTools::Math::Transformations2D;
    typedef Matrix<3,3> Matrix3x3;
    Matrix3x3 T = translation(1.0, -2.0), S = scaling(2.0, 1.0), R = rotationDegrees(45.0);
    Vector<3> p{-2,3,1};
    EXPECT_EQ(T*p, Vector<3>(-1,1,1));
    EXPECT_EQ(S*p, Vector<3>(-4,3,1));
    EXPECT_EQ(R*p, Vector<3>(-3/sqrt(2)-sqrt(2), 3/sqrt(2)-sqrt(2),1));
}

TEST(MathTest, Transformations3D)
{
    using GeometricTools::Math::Vector;
    using GeometricTools::Math::Matrix;
    using namespace GeometricTools::Math::Transformations3D;
    typedef Matrix<4,4> Matrix4x4;
    Matrix4x4 T = translation(1.0, -2.0, 3.0), S = uniformScaling(4.0), Rx = rotationXDegrees(45.0),
            Ry = rotationYDegrees(45.0), Rz = rotationZDegrees(45.0);
    Vector<4> p{-2,3,3,1};
    EXPECT_EQ(T*p, Vector<4>(-1,1,6,1));
    EXPECT_EQ(S*p, Vector<4>(-8,12,12,1));
    EXPECT_EQ(Rz*p, Vector<4>(-3/sqrt(2)-sqrt(2), 3/sqrt(2)-sqrt(2),3,1));
    EXPECT_EQ(Ry*p, Vector<4>(-3/sqrt(2)-sqrt(2), 3, 3/sqrt(2)-sqrt(2),1));
    EXPECT_EQ(Rx*p, Vector<4>(-2, 0, 3*sqrt(2),1));
}

TEST(LinearAlgebraTest, FunctionMin1D)
{
    using namespace GeometricTools::Math::NumericalOptimization;
    struct wrapper
    {
        static double f(double x)
        {
            return x*x-2*x+2;
        }

        static double f2(double x)
        {
            return x*x+2*x+2;
        }
    };
    double min, fmin;
    double a,b,c;
    double tau = 1e-10;
    std::tie(a,b,c) = getBracket(wrapper::f, 0.0, 2.0);
    std::tie(min, fmin) = goldenSearchMinimize(wrapper::f, a, b, c, tau);
    EXPECT_NEAR(min, 1.0, tau);
    EXPECT_NEAR(fmin, 1.0, tau);
    std::tie(a,b,c) = getBracket(wrapper::f2, -2.0, 0.0);
    std::tie(min, fmin) = goldenSearchMinimize(wrapper::f2, a, b, c, tau);
    EXPECT_NEAR(min, -1.0, tau);
    EXPECT_NEAR(fmin, 1.0, tau);
    std::tie(min, fmin) = goldenSearchMinimize(wrapper::f2, 0.0, 0.5, 1.0, tau);
    EXPECT_NEAR(min, 0.0, tau);
    EXPECT_NEAR(fmin, 2.0, tau);
}

TEST(ShapeTest, LinearShapesTest)
{
    using namespace GeometricTools::Primitives;
    using namespace GeometricTools::Math;
    Line<2> line({0,1}, {-1,1});
    EXPECT_EQ(line.p(), Vector<2>(0,1));
    Segment<2> seg({1,1}, {4,8});
    EXPECT_DOUBLE_EQ(seg.lengthSq(), 58.0);
    Ray<2> ray({-20,10}, {0,-1});
    EXPECT_EQ(ray.d(), Vector<2>(0,-1));
}

TEST(ShapeTest, PolylineTest)
{
    using namespace GeometricTools::Primitives;
    using namespace GeometricTools::Math;
    Polyline<3> poly;
    poly.addPoint({2,3,4});
    poly.addPoint({1,1,1});
    EXPECT_EQ(poly.vertices()[0], Vector<3>(2,3,4));
    EXPECT_EQ(poly.vertices()[1], Vector<3>(1,1,1));
    EXPECT_EQ(poly.vertices().size(), 2);
}

TEST(ShapeTest, Shapes2DTest)
{
    using namespace GeometricTools::Primitives;
    using namespace GeometricTools::Math;
    Triangle t({0,0}, {1,0}, {0,1});
    EXPECT_DOUBLE_EQ(t.area(), 0.5);
    t.addPoint({1,1});
    EXPECT_EQ(t.vertices().size(), 3);
    Rectangle r({0,0}, 3.0, 4.0);
    EXPECT_DOUBLE_EQ(r.area(), 12.0);
    r.addPoint({12,1134});
    EXPECT_EQ(r.vertices().size(), 4);
}

TEST(DistanceTest, PointToLinearTest)
{
    using namespace GeometricTools::Primitives;
    using namespace GeometricTools::Math;
    using namespace GeometricTools::Distances;
    Vector<2> p{0,1};
    Line<2> l({0,2}, {1,1});
    EXPECT_DOUBLE_EQ(Distance(l,p), 1.0/sqrt(2));
    Ray<2> r({0,2}, {1,1});
    EXPECT_DOUBLE_EQ(Distance(r,p), 1.0);
    Segment<2> s({0,2}, {0,1.5});
    EXPECT_DOUBLE_EQ(Distance(p,s), 0.5);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
