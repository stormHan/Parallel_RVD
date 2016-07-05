
#include "PolygonAction.h"

namespace P_RVD
{
	void PolygonAction::compute_weight_sum()
	{
		for (int i = 0; i < m_vertex.size(); ++i)
		{
			weight_sum += m_vertex[i].getWeight();
		}
	}

	double PolygonAction::compute_weight()
	{
		compute_weight_sum();
		return weight_sum * compute_area();
	}

	double PolygonAction::compute_area()
	{
		double area = 0.0;

		t_index v1 = 0;
		t_index v2, v3;

		for (int i = 1; i < vertex_nb - 1; ++i)
		{
			v2 = i; v3 = i + 1;
			area += compute_triangle_arae(v1, v2, v3);
		}

		return area;
	}

	Vector3d PolygonAction::compute_center()
	{
		Vector3d center;
		for (int i = 0; i < m_vertex.size(); ++i)
		{
			center += m_vertex[i].getPosition() * m_vertex[i].getWeight();
		}

		return center / weight_sum; 
	}

	double PolygonAction::compute_triangle_arae(t_index _v1, t_index _v2, t_index _v3)
	{
		Vector3d pos1 = m_vertex[_v1].getPosition();
		Vector3d pos2 = m_vertex[_v2].getPosition();
		Vector3d pos3 = m_vertex[_v3].getPosition();

		return Math::computeTriangleArea(pos1, pos2, pos3);

	}
}