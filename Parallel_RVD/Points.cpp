
#include "Points.h"

namespace P_RVD
{
	Points::Points(const Mesh& _M)
	{
		m_points = _M.meshVertices.m_Vertices;
		points_nb = m_points.size();
	}

	void Points::savePointsWithSeedStore(std::vector<SeedWeightPosition> _pos)
	{
		clear();

		points_nb = _pos.size();

		for (int i = 0; i < points_nb; ++i)
		{
			m_points.push_back(Vector3d(_pos[i].center));
		}
	}
}