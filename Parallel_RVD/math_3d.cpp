#include "math_3d.h"

Vector3f Vector3f::Cross(const Vector3f& v) const
{
	const float _x = y * v.z - z * v.y;
	const float _y = z * v.x - x * v.z;
	const float _z = x * v.y - y * v.x;

	return Vector3f(_x, _y, _z);
}

Vector3f& Vector3f::Normalize()
{
	const float Length = sqrtf(x * x + y * y + z * z);

	if (Length < 1e09 && Length > -1e09)
		return *this;

	x /= Length;
	y /= Length;
	z /= Length;

	return *this;
}

Vector3d& Vector3d::Cross(const Vector3d v) const
{
	const double _x = y * v.z - z * v.y;
	const double _y = z * v.x - x * v.z;
	const double _z = x * v.y - y * v.x;

	return Vector3d(_x, _y, _z);
}

Vector3d& Vector3d::normalize()
{
	const double Length = sqrt(x * x + y * y + z * z);

	//if (Length < 1e09 && Length > -1e09)
	//	return *this;

	x /= Length;
	y /= Length;
	z /= Length;

	return *this;
}

