#include <iostream>
#include "MatrixAlgorithms.h"
using namespace std;


int main(int argc, char *argv[])
{
	/*Vector3 test = Vector3(10.,20.,30.);
	Vector3 t = Vector3(1.,20.,50.);
	test = 2*test;
	test.Normalize();
	cout<<test[0]<<" "<<test[1]<<" "<<test[2]<<endl;
	cout<<test.Length()<<endl;*/
	/*Matrix<3,4> m = Matrix<3,4>(0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.);
	//m -= Matrix<3,4>(0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.);
	auto m2 = m.Transpose();
	auto m3 = m*m2;*/
	/*Vector3 v = Vector3(1,5,4);
	//Matrix<3,3> k = Matrix<3,3>(1.,2.,3.,4.,5.,6.,7.,8.,9.);
	Matrix<5,5> k = Matrix<5,5>(2.,3.,4.,1.,-9.,5.,4.,8.,10.,-20.,11.,29.,78.,80.,-100., 65.,-1.,0.,2.,0.,30.,-2.,-40.,7.,4.);
	k.Print();
	cout<<Determinant(k)<<endl;
	cout<<k.Norm()<<endl;*/
	Matrix<3,3> A = Matrix<3,3>(3.,-6.,-3,2.,0.,6.,-4.,7.,4.);
	Vector<3> B = Vector<3>(-3.,-22.,3.);
	Vector<3> x = solveGauss(A,B);
	x.Print();
	cout<<"------------------------------------\n\n";
	Matrix<4,4> a = Matrix<4,4>(2,1,1,0,4,3,3,1,8,7,9,5,6,7,9,8);
	Vector<4> b = Vector<4>(1.,2.,3,4);
	Vector<4> y = solveLU(a,b);
	y.Print();
	Vector<4> y2 = solveGauss(a,b);
	y2.Print();
	return 0;
}