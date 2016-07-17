
/*
	some basic math computation
*/
#include <math.h>

#include "Common.h"
#include "math_3d.h"


#ifndef H_MATH_BASICS_H
#define H_MATH_BASICS_H

namespace P_RVD
{
	
	namespace Math
	{

		/*
			compute the centor of a triangle three vertices
		*/
		Vector3d computeCenter(const Vector3d p1, const Vector3d p2, const Vector3d p3);

		/*
			compute the square of the distance of the two points
		*/
		double computeDistance(const Vector3d p1, const Vector3d p2);

		/*
			conpute the area of a triangle
			pos1, pos2, pos3 represent the position of the 3 vertex
		*/
		double computeTriangleArea(const Vector3d p1, const Vector3d p2, const Vector3d p3);
		
		/*
			Computes the centroid of a 3d triangle with weighted points.
		*/
		void computeTriangleCentroid(
			const Vector3d& p, const Vector3d& q, const Vector3d& r,
			double a, double b, double c,
			Vector3d& Vg, double& V
			);

		/*
			compute the normal with 3 vertex
		*/
		Vector3d computeNormal(const Vector3d _v1, const Vector3d _v2, const Vector3d _v3);
	}
}

#endif /* H_MATH_BASICS_H */