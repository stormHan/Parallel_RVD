
#include "PolygonAction.h"

namespace P_RVD
{
	void PolygonAction::compute_centriod()
	{
		double area = 0.0;

		t_index v1 = 0;
		t_index v2, v3;

		Vector3d pos1, pos2, pos3;
		double d1, d2, d3;
		int trianlge_nb = vertex_nb - 2;

		double total_weight = 0;
		Vector3d centriodTimesWeight;

		for (int i = 1; i < vertex_nb - 1; ++i)
		{
			v2 = i; v3 = i + 1;

			pos1 = m_vertex[v1].getPosition();
			pos2 = m_vertex[v2].getPosition();
			pos3 = m_vertex[v3].getPosition();

			d1 = m_vertex[v1].getWeight();
			d2 = m_vertex[v2].getWeight();
			d3 = m_vertex[v3].getWeight();

			Math::computeTriangleCentroid(pos1, pos2, pos3, d1, d2, d3, centriodTimesWeight, total_weight);

			current_weight += total_weight;
			current_posTimesWeight += centriodTimesWeight;

			total_weight = 0;
			centriodTimesWeight = Vector3d(0, 0, 0);
		}
		current_posTimesWeight = current_posTimesWeight / current_weight;
		
		current_weight /= trianlge_nb;
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