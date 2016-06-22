
#include "Points.h"

namespace P_RVD
{
	Points::Points(const Mesh& _M)
	{
		m_points = _M.meshVertices.m_Vertices;
		points_nb = m_points.size();
	}

}