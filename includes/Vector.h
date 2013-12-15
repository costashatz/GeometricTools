#include <iostream>
#include <limits>
#include <cmath>
#include <cstring>
using namespace std;


template<unsigned int N=2>
class Vector
{
protected:
	double values[N];
	int index;


	template<typename... Args>
	void Initialize(double h, Args&&... args)
	{
		if(index>=N)
			index=0;
		values[index] = h;
		index++;
		Initialize(forward<Args>(args)...);
	}

	template<typename... Args>
	void Initialize()
	{
		index = 0;
	}

public:	

	Vector()
	{
		memset(values,0,N*sizeof(double));
	}

	Vector(const Vector& other)
	{
		memcpy(values,other.values,N*sizeof(double));
	}

	template<typename... Args>
	Vector(double h, Args&&... args)
	{
		index = 0;
		memset(values,0,N*sizeof(double));
		Initialize(h, args...);
	}


	Vector operator+(const Vector& other)
	{
		Vector t = Vector(*this);
		for(int i=0;i<N;i++)
			t[i] += other.values[i];
		return t;
	}

	Vector operator+=(const Vector& other)
	{
		*this = *this+other;
		return *this;
	}

	Vector operator-(const Vector& other)
	{
		Vector t = Vector(*this);
		for(int i=0;i<N;i++)
			t[i] -= other.values[i];
		return t;
	}

	Vector operator-=(const Vector& other)
	{
		*this = *this-other;
		return *this;
	}

	Vector operator*(const double& other)
	{
		Vector t = Vector(*this);
		for(int i=0;i<N;i++)
			t[i] *= other;
		return t;
	}

	Vector operator*=(const double& other)
	{
		*this = *this*other;
		return *this;
	}

	Vector operator/(const double& other)
	{
		Vector t = Vector(*this);
		for(int i=0;i<N;i++)
			t[i] /= other;
		return t;
	}

	Vector operator/=(const double& other)
	{
		*this = *this/other;
		return *this;
	}

	double Length()
	{
		double s = 0.0;
		for(int i=0;i<N;i++)
			s += values[i]*values[i];
		return sqrt(s);
	}

	void Normalize()
	{
		double l = Length();
		if(l>std::numeric_limits<double>::epsilon())
			*this /= l;
	}

	double& operator()(int i) //assert legal index
	{
		return values[i];
	}

	double& operator[](int i) //assert legal index
	{
		return values[i];
	}

	double operator*(const Vector& other)
	{
		double s=0.0;
		for(int i=0;i<N;i++)
			s += values[i]*other.values[i];
		return s;
	}

	void Print()
	{
		for(int i=0;i<N;i++)
			cout<<values[i]<<" ";
		cout<<endl;
	}

	~Vector(void)
	{
	}
};


template<int N>
Vector<N> operator*(const double& a, Vector<N>& b)
{
	return b*a;
}

template<int N>
Vector<N> operator-(Vector<N>& b)
{
	return -1*b;
}

#define Vector2 Vector<2>
#define Vector3 Vector<3>
#define Vector4 Vector<4>

/*#define Vec2 Vector<2>
#define Vector2(...) Vector<2>(2,__VA_ARGS__)

#define Vec3 Vector<3>
#define Vector3(...) Vector<3>(3,__VA_ARGS__)

#define Vec4 Vector<4>
#define Vector4(...) Vector<4>(4,__VA_ARGS__)*/