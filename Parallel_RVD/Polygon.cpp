
#include "Polygon.h"

namespace P_RVD
{
	void Polygon::initialize_from_mesh_facet(
		const Mesh* _mesh, t_index _facetIdx
		)
	{
		clear();
		Facet temp_facet = _mesh->meshFacets.getFacet(_facetIdx);
		
		Vertex* v1 = add_vertex(
			Vertex(
			_mesh->meshVertices.getPoint(temp_facet.m_v1),
			1.0,
			_facetIdx
			)
		);

		Vertex* v2 = add_vertex(
			Vertex(
			_mesh->meshVertices.getPoint(temp_facet.m_v2),
			1.0,
			_facetIdx
			)
			);

		Vertex* v3 = add_vertex(
			Vertex(
			_mesh->meshVertices.getPoint(temp_facet.m_v3),
			1.0,
			_facetIdx
			)
			);
	}

	void Polygon::clip_by_plane(Polygon& _target, Points _points, t_index _i, t_index _j)
	{
		_target.clear();
		if (getVertex_nb() == 0)
		{
			return;
		}

		//get the geometic position of the i and j.
		Vector3d position_i = _points.getPoint(_i);
		Vector3d position_j = _points.getPoint(_j);

		// Compute d = n . (2m), where n is the
		// normal vector of the bisector [i, j]
		// and m the middle point of the bisector.
		double d;
		d = (position_i + position_j).cross(position_i - position_j);

		//The predecessor of the first vertex is the last vertex
		t_index prev_index_vertex = getVertex_nb() - 1;
		const Vertex* prev_vertex = &(getVertex(prev_index_vertex));
		const Vector3d prev_vertex_position = prev_vertex->getPosition();
		
		//then we compute prev_vertex_position "cross" n 
		//prev_l = prev_vertex_position . n
		double prev_l = prev_vertex_position.cross(position_i - position_j);

		P_RVD::Sign prev_status = P_RVD::geo_sgn(2.0 * prev_l - d);

		//traverse the Vertex in this Polygon
		for (t_index k = 0; k < getVertex_nb(); ++k)
		{
			const Vertex* vertex = &(getVertex(k));
			const Vector3d vertex_position = vertex->getPosition();

			double l = vertex_position.cross(position_i - position_j);
			
			P_RVD::Sign status = P_RVD::geo_sgn(2.0 * l - d);
		}

		return;
	}
}