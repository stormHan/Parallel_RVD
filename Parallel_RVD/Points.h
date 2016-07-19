
/*
	Points to store and modify the Voronoi points
*/
#include "Mesh.h"
#include "SeedStore.h"

#ifndef H_POINTS_H
#define H_POINTS_H

namespace P_RVD
{
	class Points
	{
	public:
		friend class RestrictedVoronoiDiagram;
		/*
			empty constructor
		*/
		Points(){ points_nb = 0; }
		/*
			store the Points infomation from a Mesh
		*/
		Points(const Mesh& _M);

		/*
			get the exact position of the Seed via index
		*/
		const Vector3d getPoint(t_index _t){ return m_points[_t];}

		/*
			clear the points
		*/
		void clear(){
			m_points.clear();
			points_nb = 0;
		}

		/*
			get the number of Points
		*/
		int getPointsNumber() { return points_nb; }

		/*
			set m_points with a vector from seedstore class
		*/
		void savePointsWithSeedStore(std::vector<SeedWeightPosition> _pos);

	protected:
		std::vector<Vector3d> m_points;
		int points_nb;
		
	};
}

#endif /* H_POINTS_H */