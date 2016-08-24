#include "Mesh.h"
 
namespace P_RVD
{
	/*--------------------- Line -----------------------*/
	Line::Line(t_index _begin, t_index _end)
	{
		m_begin = _begin;
		m_end = _end;
	}

	/*--------------------- Facet-----------------------*/

	Facet::Facet(t_index _v1, t_index _v2, t_index _v3)
	{
		m_v1 = _v1;
		m_v2 = _v2;
		m_v3 = _v3;
	}

	/*----------------- Mesh Vertices ------------------*/

	Vector3d MeshVertices::getPoint(t_index _nb) const
	{
		return m_Vertices[_nb];
	}
	 
	void MeshVertices::addPoint(Vector3d _newpoint)
	{
		m_nb++;
		m_Vertices.push_back(_newpoint);
	}

	MeshVertices::MeshVertices(Mesh& _mesh)
	{
		m_Vertices = _mesh.meshVertices.m_Vertices;
		m_nb = _mesh.meshVertices.m_nb;
	}

	/*----------------- Mesh Edges --------------------*/

	MeshEdges::MeshEdges(Mesh& _mesh)
	{
		m_Edges = _mesh.meshEdges.m_Edges;
	}

	/*----------------- Mesh Facets -------------------*/

	MeshFacets::MeshFacets(Mesh& _mesh)
	{
		m_Facets = _mesh.meshFacets.m_Facets;
		m_nb = _mesh.meshFacets.m_nb;
	}

	void MeshFacets::addFacet(const Facet _f)
	{
		m_nb++;
		m_Facets.push_back(_f);
	}

	void  MeshFacets::addFacet(t_index _v1, t_index _v2, t_index _v3)
	{
		m_nb++;
		m_Facets.push_back(Facet(_v1, _v2, _v3));
	}

	Facet MeshFacets::getFacet(t_index _t) const
	{
		return m_Facets[_t];
	}
	/*--------------------- Mesh ----------------------*/
	Mesh::Mesh()
		: meshVertices(*this),
		  meshEdges(*this),
		  meshFacets(*this)
	{

	}
}