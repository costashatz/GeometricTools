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

#ifndef QUADTREE_H
#define QUADTREE_H

/**
* Includes
**/
#include <Primitives/LinearShapes.h>
#include <Primitives/2D/Rectangle.h>
#include <Primitives/2D/Polygon.h>
#include <Primitives/Polyline.h>
#include <Intersections/IntersectionInfo.h>
#include <Intersections/2D/RectangleToRectangle.h>
#include <Primitives/Tools/BoundingBox.h>
#include <vector>
#include <limits>
#include <algorithm>

using std::vector;
using std::find;

namespace GeometricTools {

using Math::Vector;
using Primitives::Rectangle;
using Primitives::Polygon;
using Primitives::Polyline;
using Primitives::Segment;
using namespace Intersections;

namespace SpacePartitioning {

class QuadTree {
protected:
    Rectangle boundary_;
    QuadTree** children_;
    int level_;
    vector<Polygon> objects_;
    unsigned int max_objects_;
    unsigned int max_level_;
public:
    QuadTree(): level_(0), max_objects_(0), max_level_(0) {}

    QuadTree(const Rectangle& bounds, const int& level, const unsigned int& max_level, const unsigned int& max_objects = 1)
        : boundary_(bounds), level_(level), max_level_(max_level), max_objects_(max_objects)
    {
        children_ = new QuadTree*[4];
        for(int i=0;i<4;i++)
            children_[i] = nullptr;
    }

    bool addObject(const Polygon& obj)
    {
        if(level_ == max_level_ || !canContainObject(obj))
            return false;
        if(children_[0]!=nullptr)
        {
            for(int i=0;i<4;i++)
            {
                if(children_[i]->addObject(obj))
                    return true;
            }
            return false;
        }

        objects_.push_back(obj);
        if(objects_.size()>max_objects_)
        {
            if(!subdivide())
            {
                objects_.erase(objects_.end());
                return false;
            }
            for(int j=0;j<objects_.size();j++)
            {
                int i;
                for(i=0;i<4;i++)
                {
                    if(children_[i]->addObject(objects_[j]))
                        break;
                }
                if(i==4)
                {
                    objects_.erase(objects_.end());
                    clearChildren();
                    return false;
                }
            }
            objects_.clear();
        }
        return true;
    }

    bool removeObject(const Polygon& obj)
    {
        if(!canContainObject(obj))
            return false;
        if(children_[0]!=nullptr)
        {
            for(int i=0;i<4;i++)
            {
                if(children_[i]->removeObject(obj))
                    return true;
            }
            return false;
        }
        if(objects_.size()==0)
            return false;
        objects_.erase(find(objects_.begin(), objects_.end(), obj));
        return true;
    }

    bool queryObject(const Polygon& obj)
    {
        if(!canContainObject(obj))
            return false;
        if(children_[0]!=nullptr)
        {
            for(int i=0;i<4;i++)
            {
                if(children_[i]->queryObject(obj))
                    return true;
            }
            return false;
        }
        if(objects_.size()==0)
            return false;
        for(int i=0;i<objects_.size();i++)
        {
            Intersection2DInfo* info = intersect(Primitives::boundingBox(objects_[i]), Primitives::boundingBox(obj));
            if(info!=nullptr)
                return true;
        }
        return false;
    }

    bool queryPolyline(const Polyline<2>& poly)
    {
        Rectangle poly_bound = Primitives::boundingBox(poly);
        if(!canContainObject(poly_bound))
        if(children_[0]!=nullptr)
        {
            for(int i=0;i<4;i++)
            {
                if(children_[i]->queryPolyline(poly))
                    return true;
            }
            return false;
        }
        if(objects_.size()==0)
            return false;
        Intersection2DInfo* info;
        for(int i=0;i<objects_.size()-1;i++)
        {
            info = intersect(poly_bound, Primitives::boundingBox(objects_[i]));
            if(info != nullptr)
                return true;
        }
        return false;
    }

    Intersection2DInfo* intersectPolyline(const Polyline<2>& poly)
    {
        Rectangle poly_bound = Primitives::boundingBox(poly);
        if(!canContainObject(poly_bound))
            return nullptr;
        if(children_[0]!=nullptr)
        {
            for(int i=0;i<4;i++)
            {
                Intersection2DInfo* s = children_[i]->intersectPolyline(poly);
                if(s != nullptr)
                    return s;
            }
            return nullptr;
        }
        if(objects_.size()==0)
            return nullptr;
        Intersection2DInfo* info;
        for(int i=0;i<objects_.size()-1;i++)
        {
            info = intersect(poly_bound, Primitives::boundingBox(objects_[i]));
            if(info != nullptr)
                return info;
        }
        return nullptr;
    }

    bool full()
    {
        if(level_ == max_level_ || children_[0]==nullptr)
            return (objects_.size()==max_objects_);
        for(int i=0;i<4;i++)
        {
            if(!children_[i]->full())
                return false;
        }
        return true;
    }

    bool empty()
    {
        if(level_ == max_level_ || children_[0]==nullptr)
            return (objects_.size()==0);
        for(int i=0;i<4;i++)
        {
            if(!children_[i]->empty())
                return false;
        }
        return true;
    }

protected:
    bool subdivide()
    {
        if((level_+1)>max_level_)
            return false;
        Vector<2> p = boundary_.vertices()[0];
        Vector<2> e1 = boundary_.vertices()[1]-p;
        Vector<2> e2 = boundary_.vertices()[3]-p;
        Vector<2> center = p+(e1+e2)/2.0;
        double a = e1.length(), b = e2.length(), a_2 = a/2.0, b_2 = b/2.0;
        Vector<2> t1{-a_2, -b_2}, t2{a_2, -b_2}, t3{-a_2, b_2}, t4=-t1;

        Rectangle r = Rectangle(center+t1, a_2, b_2);
        children_[0] = new QuadTree(r, level_+1, max_level_, max_objects_);

        r = Rectangle(center+t2, a_2, b_2);
        children_[1] = new QuadTree(r, level_+1, max_level_, max_objects_);

        r = Rectangle(center+t3, a_2, b_2);
        children_[2] = new QuadTree(r, level_+1, max_level_, max_objects_);

        r = Rectangle(center+t4, a_2, b_2);
        children_[3] = new QuadTree(r, level_+1, max_level_, max_objects_);

        return true;
    }

    bool canContainObject(const Polygon& obj)
    {
        Intersection2DInfo* info = intersect(boundary_, Primitives::boundingBox(obj));
        return (info != nullptr);
    }

    void clear()
    {
        objects_.clear();
        clearChildren();
    }

    void clearChildren()
    {
        for(int i=0;i<4;i++)
        {
            if(children_[i]!=nullptr)
                children_[i]->clear();
            children_[i] = nullptr;
        }
    }
};

} }

#endif
