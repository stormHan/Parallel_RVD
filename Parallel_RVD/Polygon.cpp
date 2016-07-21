
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
		d = (position_i + position_j).dot(position_i - position_j);

		//The predecessor of the first vertex is the last vertex
		t_index prev_index_vertex = getVertex_nb() - 1;
		const Vertex* prev_vertex = &(getVertex(prev_index_vertex));
		Vector3d prev_vertex_position = prev_vertex->getPosition();
		
		//then we compute prev_vertex_position "cross" n 
		//prev_l = prev_vertex_position . n
		double prev_l = prev_vertex_position.dot(position_i - position_j);

		P_RVD::Sign prev_status = P_RVD::geo_sgn(2.0 * prev_l - d);

		//traverse the Vertex in this Polygon
		for (t_index k = 0; k < getVertex_nb(); ++k)
		{
			const Vertex* vertex = &(getVertex(k));
			Vector3d vertex_position = vertex->getPosition();

			double l = vertex_position.dot(position_i - position_j);
			
			//We compute:
			// side1(pi, pj, q) = sign(2*q.n - n.m) = sign(2*l - d)
			P_RVD::Sign status = P_RVD::geo_sgn(2.0 * l - d);

			// If status of edge extremities differ,
			// then there is an intersection.
			if (status != prev_status && (prev_status) != 0)
			{
				// create the intersection and update the Polygon
				Vertex I;
				Vector3d temp_position;

				//compute the position and weight
				double denom = 2.0 * (prev_l - l);
				double lambda1, lambda2;

				// Shit happens ! [Forrest Gump]
				if (::fabs(denom) < 1e-20) {
					lambda1 = 0.5;
					lambda2 = 0.5;
				}
				else {
					lambda1 = (d - 2.0 * l) / denom;
					// Note: lambda2 is also given
					// by (2.0*l2-d)/denom
					// (but 1.0 - lambda1 is a bit
					//  faster to compute...)
					lambda2 = 1.0 - lambda1;
				}

				temp_position.x = lambda1 * prev_vertex_position.x + lambda2 * vertex_position.x;
				temp_position.y = lambda1 * prev_vertex_position.y + lambda2 * vertex_position.y;
				temp_position.z = lambda1 * prev_vertex_position.z + lambda2 * vertex_position.z;

				//Set the position of Vertex
				I.setPosition(temp_position);

				//Set the weight of Veretex
				I.setWeight(lambda1 * prev_vertex->getWeight() + lambda2 * vertex->getWeight());

				//???? »¹Ã»Åª¶®
				if (status > 0)
				{
					I.copy_edge_from(*prev_vertex);
					I.setSeed(signed_t_index(_j));
				}
				else
				{
					I.setEdgeType(Intersection);
					I.setSeed(vertex->getSeed());
				}
				_target.add_vertex(I);
			}

			if (status > 0)
			{
				_target.add_vertex(*vertex);
			}
			prev_vertex = vertex;
			prev_vertex_position = vertex_position;
			prev_status = status;
			prev_l = l;
			prev_index_vertex = k;
		}

		return;
	}
}