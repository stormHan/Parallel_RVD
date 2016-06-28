//---------------------------------------
//
//    stormhan  math_3d.h
//
//--------------------------------------


#ifndef H_MATH_3D_H					//防止重复定义
#define H_MATH_3D_H

#include <iostream>
#ifdef WIN32						//为了增强在不同平台上的兼容性
#define _USE_MATH_DEFINES
#include <cmath>

#else
#include <math.h>
#endif

#ifndef M_PI
#define M_PI 3.1415926
#endif

#define ToRadian(x) (float)((x) * M_PI / 180.0f)   
#define ToDegree(x) (float)((x) * 180.0f / M_PI)


struct Vector2i
{
	int x;
	int y;
};

struct Vector3i
{
	int x;
	int y;
	int z;

	Vector3i()
	{

	}

	Vector3i(int _x, int _y, int _z)
	{
		x = _x;
		y = _y;
		z = _z;
	}

	Vector3i(int f)
	{
		x = y = z = f;
	}

	Vector3i& operator+=(const Vector3i& r)
	{
		x += r.x;
		y += r.y;
		z += r.z;

		return *this;
	}

	Vector3i& operator-=(const Vector3i& r)
	{
		x -= r.x;
		y -= r.y;
		z -= r.z;

		return *this;
	}

	Vector3i& operator*=(const Vector3i& r)
	{
		x *= r.x;
		y *= r.y;
		z *= r.z;

		return *this;
	}
};

struct Vector2f
{
	float x;
	float y;

	Vector2f()
	{

	}

	Vector2f(float _x, float _y)
	{
		x = _x;
		y = _y;
	}

};

struct Vector3f
{
	float x;
	float y;
	float z;

	Vector3f()
	{

	}

	Vector3f(float _x, float _y, float _z)
	{
		x = _x;
		y = _y;
		z = _z;
	}

	Vector3f(float f)
	{
		x = y = z = f;
	}

	Vector3f& operator+=(const Vector3f& r)
	{
		x += r.x;
		y += r.y;
		z += r.z;

		return *this;
	}

	Vector3f& operator-=(const Vector3f& r)
	{
		x -= r.x;
		y -= r.y;
		z -= r.z;

		return *this;
	}

	Vector3f& operator*=(const Vector3f& r)
	{
		x *= r.x;
		y *= r.y;
		z *= r.z;

		return *this;
	}

	Vector3f& Normalize();

	Vector3f Cross(const Vector3f& v) const;

	//void Rotate(float Angle, const Vector3f& Axis);

	void Print() const
	{
		printf("(%.02f,  %.02f, %.02f)", x, y, z);
	}

	void Rotate(float Angle, const Vector3f& Axe);
};

struct Vector4f
{
	float x;
	float y;
	float z;
	float w;

	Vector4f()
	{
	}

	Vector4f(float _x, float _y, float _z, float _w)
	{
		x = _x;
		y = _y;
		z = _z;
		w = _w;
	}

	Vector4f(float f)
	{
		x = y = z = w = f;
	}

	void Print() const										//C++中函数后面的const 表示该函数不影响成员变量的值
	{
		printf("(%.02f, %.02f, %.02f, %.02f)", x, y, z, w);
	}

	Vector3f to3f() const
	{
		Vector3f v(x, y, z);
		return v;
	}
};

struct Vector3d
{
	double x;
	double y;
	double z;

	Vector3d()
	{

	}

	Vector3d(double _x, double _y, double _z)
	{
		x = _x;
		y = _y;
		z = _z;
	}

	Vector3d(double f)
	{
		x = y = z = f;
	}

	Vector3d& operator+=(Vector3d& r)
	{
		x += r.x;
		y += r.y;
		z += r.z;

		return *this;
	}

	Vector3d& operator-=(const Vector3d& r)
	{
		x -= r.x;
		y -= r.y;
		z -= r.z;

		return *this;
	}

	Vector3d& operator =(const Vector3d& r)
	{
		x = r.x;
		y = r.y;
		z = r.z;

		return *this;
	}

	double cross(const Vector3d r) const
	{
		return x * r.x + y * r.y + z * r.z;
	}
};



inline Vector3d operator+(const Vector3d& l, const Vector3d& r)
{
	Vector3d Ret(l.x + r.x,
				 l.y + r.y,
				 l.z + r.z);

	return Ret;
}

inline Vector3d operator-(const Vector3d& l, const Vector3d& r)
{
	Vector3d Ret(l.x - r.x,
				 l.y - r.y,
				 l.z - r.z);
	return Ret;
}

inline Vector3d operator*(const Vector3d& l, const float f)
{
	Vector3d Ret(l.x * f,
				 l.y * f,
				 l.z * f);
	return Ret;
}

#endif /* H_MATH_3D_H*/