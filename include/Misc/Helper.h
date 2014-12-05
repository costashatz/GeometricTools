#ifndef HELPER_H
#define HELPER_H

namespace LinearAlgebraTools { namespace Helper {

const double Pi = 3.14159265358979323846264338327950288419716939937510;
const double   TwoPi     = Pi * 2;
const double   HalfPi    = Pi / 2;
const double   QuarterPi = Pi / 4;

double radiansToDegrees(double rad)
{
    return rad*(180.0/Pi);
}

double degreesToRadians(double deg)
{
    return deg*(Pi/180.0);
}

} }

#endif
