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

#include <iostream>
#include <geometric_tools/Primitives/2D/PlaneCurve.h>
#include <geometric_tools/Primitives/1D/CardinalPiecewiseCurve.h>
#include <geometric_tools/Primitives/1D/HermitePiecewiseCurve.h>
#include <geometric_tools/Primitives/Tools/CurveSubdivision.h>
using namespace std;

using namespace GeometricTools::Primitives;

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
    Polyline<2> poly = uniformSubdivision(hermite, 10);
    return 0;
}
