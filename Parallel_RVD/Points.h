
/*
	Points to store and modify the Voronoi points
*/
#include "Mesh.h"


#ifndef H_POINTS_H
#define H_POINTS_H

namespace P_RVD
{
	class Points
	{
	public:
		friend class RestrictedVoronoiDiagram;
		/*
			store the Points infomation from a Mesh
		*/
		Points(const Mesh& _M);

		/*
			get the exact position of the Seed via index
		*/
		const Vector3d getPoint(t_index _t){ return m_points[_t];}

	protected:
		std::vector<Vector3d> m_points;
		int points_nb;
		
	};
}

#endif /* H_POINTS_H */