#include <iostream>
#include <Primitives/2D/PlaneCurve.h>
#include <Primitives/1D/CardinalPiecewiseCurve.h>
#include <Primitives/1D/HermitePiecewiseCurve.h>
#include <Primitives/Tools/CurveSubdivision.h>
using namespace std;

using namespace LinearAlgebraTools::Primitives;

//SIMPLE EXAMPLE

int main(int argc, char *argv[])
{
    // Create Catmull-Rod piecewise 2D curve
    PlaneCurve<CardinalPiecewiseCurve,CardinalPiecewiseCurve> catmull_rod(0.0);
    catmull_rod.addPoint({1,2});
    catmull_rod.addPoint({2,7});
    catmull_rod.addPoint({3,0});
    catmull_rod.addPoint({2,-3});
    catmull_rod.addPoint({0,-7});
    cout<<"CATMULL-ROD:\n\n";
    // Get coefficients for each axis
    for(int i=0;i<2;i++)
    {
        auto coef = catmull_rod.coeff(i);
        char x;
        if(i==0)
            x = 'x';
        else
            x = 'y';
        for(int j=0;j<coef.size();j++)
        {
            int ind = (j==1)? j-1 : j;
            if((j%4)!=0 && coef[j]>0.0 && coef[ind]>0.0)
                cout<<"+ ";
            if((j%4)==0)
                cout<<x<<j/4<<": ";
            if(coef[j] != 0.0)
                cout<<coef[j];
            if(coef[j] != 0.0 && (3-(j%4))!=0.0)
                cout<<"*u^"<<3-(j%4)<<" ";
            if(j%4==3)
                cout<<endl;
        }
    }
    // Here you could use gnu-Plot or Boost-GIL to plot the curve

    // Create Hermite piecewise 2D curve
    PlaneCurve<HermitePiecewiseCurve,HermitePiecewiseCurve> hermite;
    hermite.addDotPoint({1,2,0,1});
    hermite.addDotPoint({2,7,1,2});
    hermite.addDotPoint({3,0,-1,-2});
    hermite.addDotPoint({2,-3,10,-8});
    hermite.addDotPoint({0,-7,3,0});
    hermite.addDotPointInPlace({1.5, 2.5, 0.1, 0.5}, 1);
    cout<<"\n\nHERMITE:\n\n";
    // Get coefficients for each axis
    for(int i=0;i<2;i++)
    {
        auto coef = hermite.coeff(i);
        char x;
        if(i==0)
            x = 'x';
        else
            x = 'y';
        for(int j=0;j<coef.size();j++)
        {
            int ind = (j==1)? j-1 : j;
            if((j%4)!=0 && coef[j]>0.0 && coef[ind]>0.0)
                cout<<"+ ";
            if((j%4)==0)
                cout<<x<<j/4<<": ";
            if(coef[j] != 0.0)
                cout<<coef[j];
            if(coef[j] != 0.0 && (3-(j%4))!=0.0)
                cout<<"*u^"<<3-(j%4)<<" ";
            if(j%4==3)
                cout<<endl;
        }
    }
    // Subdivide curve to polyline using N segments
    Polyline poly = uniformSubdivision(hermite, 10);
    return 0;
}
