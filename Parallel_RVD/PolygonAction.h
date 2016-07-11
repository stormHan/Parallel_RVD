/*
	a class used to compute the polygon infomation.	
*/

#ifndef H_POLYGONACTION_H
#define H_POLYGONACTION_H

#include "Polygon.h"
#include "math_3d.h"
#include "Math_basics.h"

namespace P_RVD
{
	class PolygonAction
	{

	public:

		PolygonAction(Polygon* _p)
		{
			m_polygon = _p;
			weight_sum = 0.0;
		}

		PolygonAction(){
			m_polygon = NULL;
			weight_sum = 0.0;
		}

		/*
			get the computing Polygon
		*/
		Polygon* getPolygon()
		{
			return m_polygon;
		}

		/*
			clear the PolygonAction to get prepared for the next action.
		*/
		void clear()
		{
			m_polygon = NULL;
			weight_sum = 0.0; 
			m_vertex.clear();
			vertex_nb = m_vertex.size();

			current_weight = 0.0;
			current_posTimesWeight = Vector3d(0.0, 0.0, 0.0);
		}

		/*
			set the computing
		*/
		void setPolygon(Polygon* _p)
		{
			m_polygon = _p;

			vertex_nb = m_polygon->getVertex_nb();

			for (int i = 0; i < vertex_nb; ++i)
			{
				m_vertex.push_back(m_polygon->getVertex(i));
			}
		}

		/*
			compute the total weight and centriod times weight
		*/
		void compute_centriod();

		/*
			compute the weight of the m_polygon
		*/
		double compute_weight();

		/*
			compute the area of the m_polygin
		*/
		double compute_area();

		/*
			compute the center of the m_polygon
		*/
		Vector3d compute_center();

		/*
			compute a triangle's area by index
			_v1, _v2, _v3 represent the indices of m_vertex
		*/
		double compute_triangle_arae(t_index _v1, t_index _v2, t_index _v3);

		/*
			get the current_position
		*/
		Vector3d getCurrentPosition(){ return current_posTimesWeight; }

		/*
			get the current weight
		*/
		double getCurrentWeight(){ return current_weight; }

	protected:
		Polygon* m_polygon;
		std::vector<Vertex> m_vertex;
		int vertex_nb;

		Vector3d current_posTimesWeight;
		double current_weight;

		double weight_sum;
	};

}

#endif /* H_POLYGONACTION_H */