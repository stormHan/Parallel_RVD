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
