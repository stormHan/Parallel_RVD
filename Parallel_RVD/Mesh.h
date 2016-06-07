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
	class Line
	{
	public:
		Line(t_index _begin, t_index _end);

	private:
		t_index m_begin;
		t_index m_end;
	};

	class Facet
	{
	public:
		Facet(t_index _v1, t_index _v2, t_index _v3);

	private:
		t_index m_v1;
		t_index m_v2;
		t_index m_v3;
	};

	class Mesh;
	/*----------------- Mesh Vertices ------------------*/
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

	private:
		std::vector<Line> m_Edges;
	};

	class MeshFacets
	{
	public:
		MeshFacets(Mesh& _mesh);

		/*
			add a Facet to m_Facets
		*/
		void addFacet(const Facet _f);

		/*
			add a Facet via three indice of m_Vertices
		*/
		void addFacet(t_index _v1, t_index _v2, t_index _v3);
		
	private:
		std::vector<Facet> m_Facets;
		unsigned int m_nb;
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