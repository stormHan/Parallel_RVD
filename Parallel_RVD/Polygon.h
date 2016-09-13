/*
	a new data stucture to store the information
	about the computation of RVD.
	we compute the intersection of a polygon and a seed(cell)
*/

#ifndef H_POLYGON_H
#define H_POLYGON_H

#include <vector>

#include "RVD_Vertex.h"
#include "Common.h"
#include "Mesh.h"
#include "Points.h"

namespace P_RVD
{
	class Vertex;

	class Polygon
	{
	public:
		Polygon(){}

		Polygon(Polygon& _p)
		{
			m_vertex = _p.m_vertex;
		}
		/*
			get the number of vertex in this Polygon
		*/
		t_index getVertex_nb() const
		{
			return (t_index)m_vertex.size();
		}
		
		/*
			get a Vertex in m_vertex by indice.
		*/
		const Vertex& getVertex(t_index _t) const
		{
			if (_t >= getVertex_nb())
			{
				printf("wrong index to get the vertex : %d", _t);
				exit(0);
			}
			return m_vertex[_t];
		}

		/*
			get the index of the next Vertex of this Vertex
			if this Vertex is the last one of this Polygon
			return 0 ( the first one )
		*/
		t_index getNextVertexIndex(t_index _t)
		{
			if (_t >= getVertex_nb())
			{
				printf("wrong index to get the next vertex index : %d", _t);
				exit(0);
			}
			return (_t == getVertex_nb() - 1) ? 0 : (_t + 1);
		}

		/*
			Add a Vertex to this Polygon.
			return the address of the stored Vertex.
		*/
		Vertex* add_vertex(const Vertex& v) {
			m_vertex.push_back(v);
			return &*(m_vertex.rbegin());
		}

		/*
			clear the Polygon to empty
		*/
		void clear()
		{
			m_vertex.resize(0);
		}

		/*
			clear the Polygon to a appointed size
		*/
		void clear(t_index _size)
		{
			m_vertex.resize(_size);
		}

		void initialize_from_mesh_facet(
			const Mesh* _mesh, t_index _facetIdx
			);

		/*
			be clipped by a bisector
			points : the seeds set.
			i : the core seed
			j : the seed around the core seed i.
		*/
		void clip_by_plane(Polygon& _target, Points _points, t_index _i, t_index _j);

		/*
			print the polygon
		*/
		void show_polygon();

	protected:
		std::vector<Vertex> m_vertex;
	};
}
#endif /* H_POLYGON_H */