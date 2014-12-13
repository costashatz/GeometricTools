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
    EXPECT_EQ(c.normalized().length(), 1.0);
    d.normalize();
    EXPECT_EQ(d.length(), 1.0);
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
    EXPECT_EQ(A.norm(), sqrt(3));
}

TEST(MathTest, Transformations2D)
{
    using GeometricTools::Math::Vector;
    using GeometricTools::Math::Matrix;
    using namespace GeometricTools::Math::Transformations2D;
    typedef Matrix<3,3> Matrix3x3;
    Matrix3x3 T = translation(1.0, -2.0), S = scaling(2.0, 1.0), R = rotationDegrees(45.0);
    Vector<3> p{-2,3,1.0};
    EXPECT_EQ(T*p, Vector<3>(-1,1,1));
    EXPECT_EQ(S*p, Vector<3>(-4,3,1));
    EXPECT_EQ(R*p, Vector<3>(-3/sqrt(2)-sqrt(2), 3/sqrt(2)-sqrt(2),1));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
