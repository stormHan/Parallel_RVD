#include "Mesh.h"
 
namespace P_RVD
{
	/*----------------- Mesh Vertices ------------------*/

	Vector3f MeshVertices::getPoint(t_index _nb)
	{
		return m_Vertices[_nb];
	}

	void MeshVertices::addPoint(Vector3f _newpoint)
	{
		m_nb++;
		m_Vertices.push_back(_newpoint);
	}

	MeshVertices::MeshVertices(Mesh& _mesh)
	{
		m_Vertices = _mesh.meshVertices.m_Vertices;
	}

	/*----------------- Mesh Edges --------------------*/

	MeshEdges::MeshEdges(Mesh& _mesh)
	{
		
	}

	/*----------------- Mesh Edges --------------------*/

	MeshFacets::MeshFacets(Mesh& _mesh)
	{

	}
	/*--------------------- Mesh ----------------------*/
	Mesh::Mesh()
		: meshVertices(*this),
		  meshEdges(*this),
		  meshFacets(*this)
	{

	}
}