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



struct Quaternion
{
	float x, y, z, w;

	Quaternion(float _x, float _y, float _z, float _w);

	void Normalize();

	Quaternion Conjugate();

	Vector3f ToDegrees();
};

class Matrix4f
{
public:
	float m[4][4];

	Matrix4f()
	{
	}

	Matrix4f(float a00, float a01, float a02, float a03,
		float a10, float a11, float a12, float a13,
		float a20, float a21, float a22, float a23,
		float a30, float a31, float a32, float a33)
	{
		m[0][0] = a00; m[0][1] = a01; m[0][2] = a02; m[0][3] = a03;
		m[1][0] = a10; m[1][1] = a11; m[1][2] = a12; m[1][3] = a13;
		m[2][0] = a20; m[2][1] = a21; m[2][2] = a22; m[2][3] = a23;
		m[3][0] = a30; m[3][1] = a31; m[3][2] = a32; m[3][3] = a33;
	}


	Matrix4f Transpose() const
	{
		Matrix4f n;

		for (unsigned int i = 0; i < 4; i++) {
			for (unsigned int j = 0; j < 4; j++) {
				n.m[i][j] = m[j][i];
			}
		}

		return n;
	}


	inline void InitIdentity()
	{
		m[0][0] = 1.0f; m[0][1] = 0.0f; m[0][2] = 0.0f; m[0][3] = 0.0f;
		m[1][0] = 0.0f; m[1][1] = 1.0f; m[1][2] = 0.0f; m[1][3] = 0.0f;
		m[2][0] = 0.0f; m[2][1] = 0.0f; m[2][2] = 1.0f; m[2][3] = 0.0f;
		m[3][0] = 0.0f; m[3][1] = 0.0f; m[3][2] = 0.0f; m[3][3] = 1.0f;
	}

	inline Matrix4f operator*(const Matrix4f& Right) const
	{
		Matrix4f Ret;

		for (unsigned int i = 0; i < 4; i++) {
			for (unsigned int j = 0; j < 4; j++) {
				Ret.m[i][j] = m[i][0] * Right.m[0][j] +
					m[i][1] * Right.m[1][j] +
					m[i][2] * Right.m[2][j] +
					m[i][3] * Right.m[3][j];
			}
		}

		return Ret;
	}

	Vector4f operator*(const Vector4f& v) const
	{
		Vector4f r;

		r.x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3] * v.w;
		r.y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3] * v.w;
		r.z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3] * v.w;
		r.w = m[3][0] * v.x + m[3][1] * v.y + m[3][2] * v.z + m[3][3] * v.w;

		return r;
	}

	operator const float*() const
	{
		return &(m[0][0]);
	}

	void Print() const
	{
		for (int i = 0; i < 4; i++) {
			printf("%f %f %f %f\n", m[i][0], m[i][1], m[i][2], m[i][3]);
		}
	}


};

inline Vector3f operator+(const Vector3f& l, const Vector3f& r)
{
	Vector3f Ret(l.x + r.x,
				 l.y + r.y,
				 l.z + r.z);

	return Ret;
}

inline Vector3f operator-(const Vector3f& l, const Vector3f& r)
{
	Vector3f Ret(l.x - r.x,
				 l.y - r.y,
				 l.z - r.z);
	return Ret;
}

inline Vector3f operator*(const Vector3f& l, const float f)
{
	Vector3f Ret(l.x * f,
				 l.y * f,
				 l.z * f);
	return Ret;
}

inline Vector4f operator/(const Vector4f& l, const float f)
{
	if (f < 1e09 && f > -1e09)
	{
		Vector4f ret(0.0);
		return ret;
	}

	Vector4f Ret(l.x / f,
				 l.y / f,
				 l.z / f,
				 l.w / f);
	return Ret;
}

Quaternion operator*(const Quaternion& l, const Quaternion& r);

Quaternion operator*(const Quaternion& q, const Vector3f& v);
#endif /* H_MATH_3D_H*/