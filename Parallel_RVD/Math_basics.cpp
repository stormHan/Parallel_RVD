
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

		float computeDistance(const Vector3d p1, const Vector3d p2)
		{
			return (p1.x - p2.x) * (p1.x - p2.x)
				+ (p1.y - p2.y) * (p1.y - p2.y)
				+ (p1.z - p2.z) * (p1.z - p2.z);
		}
	}
}