/*
	RVD.h
	used to Compute Restricted Voronori Diagram
*/


#ifndef H_RVD_H
#define H_RVD_H
#include <algorithm>
#include <map>

#include "Common.h"
#include "Mesh.h"
#include "Points.h"
#include "Polygon.h"
#include "PolygonAction.h"

#include "SeedStore.h"

#include "Math_basics.h"

namespace P_RVD
{
	class RestrictedVoronoiDiagram
	{
	public:

		/*
			constructor function of this class
			since the computation of the RVD is base on Voronoi cells and Facets
			this class is bound with 2 object (Mesh, Points)
		*/
		RestrictedVoronoiDiagram(Mesh* _M, Points* _P);

		/*
			compute the RVD 
		*/
		bool compute_RVD();

		/*
			get the center of a triangle
		*/
		Vector3d computeCenter(const Vector3i index_triangle);
		Vector3d computeCenter(const t_index t1, const t_index t2, const t_index t3);

		/*
			find the nearest n points in p_Points of an index of the Seeds
		*/

		std::vector<int> findNearestPoints(const t_index _t, int _n);

		/*
			find the nearest n points in p_Points of a point
		*/
		std::vector<int> findNearestPoints(const Vector3d _center, int _n);

		/*
			Compute the intersection between the Voronoi Cell of a  seed and a facet.
		*/
		Polygon* intersect_cell_facet(t_index seed, Polygon& F);

		/*
			a polygon clipped by a plane.
			a plane (bisector) is defined by 2 point i and j.

		*/
		void clip_by_plane(Polygon& ping, Polygon& pong, t_index i, t_index j);

		/*
			Swap two pointers between two polygons
		*/
		void swap_polygons(Polygon*& ping, Polygon*& pong);

		/*
			deconstruction
		*/
		~RestrictedVoronoiDiagram()
		{
			delete polygonHandler;
		}


	private:
		Mesh* p_Mesh;
		Points* p_Points;

		Vector3d current_center;
		std::vector<int> current_near_points;

		t_index face_begin;
		t_index face_end;

		unsigned int seeds_n;
		t_index current_seed;

		Polygon* current_polygon;
		/*
			as a buffer the store the result after clipping
		*/
		Polygon polygon_buffer;
		/*
			empty polygon
		*/
		Polygon empty_polygon;
		/*
			the neighbors betweent a seed
		*/
		std::vector<int> seed_neighbors;

		/*
			handle the clipped polygon information
		*/
		PolygonAction* polygonHandler;

		/*
			temp weight and center of the clipped polygon
		*/
		double polygon_weight;
		Vector3d polygon_center;

		/*
			used for updating the information of seeds.
		*/
		SeedStore seedsUpdater;

		std::vector<int> polygon_nb;

	};
}

#endif /* H_RVD_H */