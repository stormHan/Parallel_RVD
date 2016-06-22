/*
	RVD.h
	used to Compute Restricted Voronori Diagram
*/


#ifndef H_RVD_H
#define H_RVD_H
#include <algorithm>
#include <map>

#include "Mesh.h"
#include "Points.h"

#include "Math_basics.h"

namespace P_RVD
{
	class RestrictedVoronoiDiagram
	{
	public:

		/*
			constructor function of this class
			since the computation of the RVD is base on Voronoi cells and Facets
			this class is bound with 2 object (Mesh, Points)
		*/
		RestrictedVoronoiDiagram(Mesh* _M, Points* _P);

		/*
			compute the RVD 
		*/
		bool compute_RVD();

		/*
			get the center of a triangle
		*/
		Vector3d computeCenter(const Vector3i index_triangle);
		Vector3d computeCenter(const t_index t1, const t_index t2, const t_index t3);

		/*
			find the nearest n points in p_Points of a point
		*/
		std::vector<int> findNearestPoints(const Vector3d _center, int _n);

	private:
		Mesh* p_Mesh;
		Points* p_Points;

		Vector3d current_center;
		std::vector<int> current_near_points;
	};
}

#endif /* H_RVD_H */