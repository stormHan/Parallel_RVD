
#include "Math_basics.h"

namespace P_RVD
{
	namespace Math
	{
		Vector3d computeCenter(const Vector3d p1, const Vector3d p2, const Vector3d p3)
		{
			return Vector3d((p1.x + p2.x + p3.x) / 3,
				(p1.y + p2.y + p3.y) / 3,
				(p1.z + p2.z + p3.z) / 3);
		}

		double computeDistance(const Vector3d p1, const Vector3d p2)
		{
			return sqrt((p1.x - p2.x) * (p1.x - p2.x)
				+ (p1.y - p2.y) * (p1.y - p2.y)
				+ (p1.z - p2.z) * (p1.z - p2.z));
		}

		double computeTriangleArea(const Vector3d p1, const Vector3d p2, const Vector3d p3)
		{
			double a = computeDistance(p1, p2);
			double b = computeDistance(p2, p3);
			double c = computeDistance(p1, p3);

			//Heron's Formula to compute the area of the triangle
			double p = (a + b + c) / 2;
			
			if (p >= a && p >= b && p >= c)
				return sqrt(p * (p - a) * (p - b) * (p - c));
			else
				fprintf(stderr, "the three point cannot construct a triangle!\n");
			return 0.0;
		}
	}
}