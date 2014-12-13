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

#ifndef PIECEWISE_CURVE_H
#define PIECEWISE_CURVE_H

/**
* Includes
**/
#include <Primitives/1D/Curve.h>

#include <Math/Vector.h>
#include <Math/Matrix.h>
#include <vector>

using std::vector;


namespace GeometricTools {

using Math::Matrix;
using Math::Vector;

namespace Primitives {

/**
* PiecewiseCurve Class
* Abstract 1D piecewise curve - base class for all 1D piecewise curves
**/
class PiecewiseCurve : public Curve
{
protected:
    vector<Curve*> curves;
public:
    PiecewiseCurve() {}

    virtual vector<double> coeff()
    {
        if(curves.size()==0)
            return vector<double>();
        vector<double> coef = curves[0]->coeff();
        for(int i=1;i<curves.size();i++)
        {
            auto tmp = curves[i]->coeff();
            coef.insert(coef.end(), tmp.begin(), tmp.end());
        }
        return coef;
    }

    double getPoint(const double &u)
    {
        // Assumes that all the piecewise curve is in u[0,1]
        int c = getCurveFromU(u*curves.size());
        if(c==-1)
            return 0.0;
        return curves[c]->getPoint(u*curves.size()-c);
    }
protected:
    int getCurveFromU(const double& u)
    {
        if(curves.size()==0)
            return -1;
        int c = int(u);
        if(c>=curves.size())
            c = curves.size()-1;
        return c;
    }
};

} }


#endif
