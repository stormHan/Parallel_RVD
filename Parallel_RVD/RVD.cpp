
#include "RVD.h"

namespace P_RVD
{
	RestrictedVoronoiDiagram::RestrictedVoronoiDiagram(Mesh* _M, Points* _P)
	{
		p_Mesh = _M;
		p_Points = _P;

		seeds_n = 5;
	}

	Vector3d RestrictedVoronoiDiagram::computeCenter(const Vector3i index_triangle)
	{
		Vector3d p1 = p_Mesh->meshVertices.getPoint(index_triangle.x - 1);
		Vector3d p2 = p_Mesh->meshVertices.getPoint(index_triangle.y - 1);
		Vector3d p3 = p_Mesh->meshVertices.getPoint(index_triangle.z - 1);

		return Math::computeCenter(p1, p2, p3);
	}

	Vector3d RestrictedVoronoiDiagram::computeCenter(const t_index t1, const t_index t2, const t_index t3)
	{
		return computeCenter(Vector3i((int)t1, (int)t2, (int)t3));
	}

	std::vector<int> RestrictedVoronoiDiagram::findNearestPoints(const Vector3d _center, int _n)
	{
		std::map<double, t_index> dis_index_map;
		std::vector<int> nearest_N_index;
		std::vector<double> temp_dis;
		double t;
		for (int i = 0; i < p_Points->points_nb; ++i)
		{
			t = Math::computeDistance(p_Points->m_points[i], _center);
			temp_dis.push_back(t);

			while (dis_index_map.count(t) != 0)
				t += 0.00000000001;
			dis_index_map.insert(std::pair<double, t_index>(t, i));
		}
		
		std::sort(temp_dis.begin(), temp_dis.end());

		for (int i = 0; i < _n; ++i)
		{
			nearest_N_index.push_back(dis_index_map[temp_dis[i]]);
		}
		return nearest_N_index;
	}

	bool RestrictedVoronoiDiagram::compute_RVD()
	{
		face_begin = 0;
		face_end = p_Mesh->meshFacets.getFacetsNumber();

		/*
			compute the RVD of each facet and each seed
		*/
		for (t_index t = face_begin; t < face_end; ++t)
		{
			Facet temp_facet = p_Mesh->meshFacets.getFacet(t);
			current_center = computeCenter(temp_facet.m_v1, temp_facet.m_v2, temp_facet.m_v3);

			current_near_points = findNearestPoints(current_center, seeds_n);

			for (t_index i = 0; i < seeds_n; ++i)
			{
				current_seed = current_near_points[i];
			}
		}
		return true;
	}
}