#ifndef HELPER_H
#define HELPER_H
#include <iostream>
#include <limits>
#define minimum(x,y) (x>y)?x:y
using namespace std;

const int     MaxInt    = (std::numeric_limits<int>::max)();
const double  MaxDouble = (std::numeric_limits<double>::max)();
const double  MinDouble = (std::numeric_limits<double>::min)();
const float   MaxFloat  = (std::numeric_limits<float>::max)();
const float   MinFloat  = (std::numeric_limits<float>::min)();

const double Pi = 3.14159265358979323846264338327950288419716939937510;
const double   TwoPi     = Pi * 2;
const double   HalfPi    = Pi / 2;
const double   QuarterPi = Pi / 4;

double RadiansToDegrees(double rad)
{
	return rad*(180.0/Pi);
}

double DegreesToRadians(double deg)
{
	return deg*(Pi/180.0);
}

#endif