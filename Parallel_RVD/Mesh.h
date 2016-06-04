/*
	***Data Structure Mesh***

	Store the infomation of 
		Input	Model;
		Input	Points;
		Output	Model;

	Include:
		Mesh Vertex;
		Mesh Edges;
		Mesh Facets;

*/

#ifndef H_MESH_H
#define H_MESH_H

#include "math_3d.h"
#include "Common.h"

#include <vector>

namespace P_RVD
{
	class Mesh;

	class MeshVertices
	{
	public:
		MeshVertices(Mesh& _mesh);
		

		/*
			get Point Position via indice
		*/
		
		Vector3f getPoint(t_index _nb);
		
		/*
			add a Point to m_Vertices
		*/
		void addPoint(Vector3f _newpoint);

	private:
		std::vector<Vector3f> m_Vertices;
		unsigned int m_nb = 0;
	};

	class MeshEdges
	{
	public:
		MeshEdges(Mesh& _mesh);
	};

	class MeshFacets
	{
	public:
		MeshFacets(Mesh& _mesh);
		
	private:

	};



	class Mesh
	{
	public:
		MeshVertices	meshVertices;
		MeshEdges		meshEdges;
		MeshFacets		meshFacets;

		
		Mesh();
	private:

	};

}
#endif /* H_MESH_H */