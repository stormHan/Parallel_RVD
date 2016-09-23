
#include "RVD.h"

namespace P_RVD
{
	RestrictedVoronoiDiagram::RestrictedVoronoiDiagram(Mesh* _M, Points* _P)
	{
		p_Mesh = _M;
		p_Points = _P;

		seeds_n = 20;
		polygonHandler = new PolygonAction();

		seedsUpdater.setSeedsNumber(p_Points->points_nb);
		seedsUpdater.setPositionVector();

		polygon_nb = std::vector<int>(p_Points->points_nb, 0);
	}

	Vector3d RestrictedVoronoiDiagram::computeCenter(const Vector3i index_triangle)
	{
		Vector3d p1 = p_Mesh->meshVertices.getPoint(index_triangle.x);
		Vector3d p2 = p_Mesh->meshVertices.getPoint(index_triangle.y);
		Vector3d p3 = p_Mesh->meshVertices.getPoint(index_triangle.z);

		return Math::computeCenter(p1, p2, p3);
	}

	Vector3d RestrictedVoronoiDiagram::computeCenter(const t_index t1, const t_index t2, const t_index t3)
	{
		return computeCenter(Vector3i((int)t1, (int)t2, (int)t3));
	}

	std::vector<int> RestrictedVoronoiDiagram::findNearestPoints(const t_index _t, int _n)
	{
		std::map<double, t_index> dis_index_map;
		std::vector<int> nearest_N_index;
		std::vector<double> temp_dis;

		Vector3d _center = p_Points->getPoint(_t);

		double t;
		for (int i = 0; i < p_Points->points_nb; ++i)
		{
			t = Math::computeDistance(p_Points->m_points[i], _center);

			while (dis_index_map.count(t) != 0)
				t += 0.00000000001;
			dis_index_map.insert(std::pair<double, t_index>(t, i));

			temp_dis.push_back(t);
		}

		std::sort(temp_dis.begin(), temp_dis.end());

		for (int i = 0; nearest_N_index.size() < _n; ++i)
		{
			int tmp = dis_index_map[temp_dis[i]];
			if (tmp != _t)
				nearest_N_index.push_back(tmp);
		}
		return nearest_N_index;
	}

	std::vector<int> RestrictedVoronoiDiagram::findNearestPoints(const Vector3d _center, int _n)
	{
		std::map<double, t_index> dis_index_map;
		std::vector<int> nearest_N_index;
		std::vector<double> temp_dis;

		double t;
		for (int i = 0; i < p_Points->points_nb; ++i)
		{
			t = Math::computeDistance(p_Points->m_points[i], _center);
			 
			while (dis_index_map.count(t) != 0)
				t += 0.00000000001;
			dis_index_map.insert(std::pair<double, t_index>(t, i));

			temp_dis.push_back(t);
		}
		
		std::sort(temp_dis.begin(), temp_dis.end());

		for (int i = 0; i < _n; ++i)
		{
			nearest_N_index.push_back(dis_index_map[temp_dis[i]]);
		}
		return nearest_N_index;
	}

	bool RestrictedVoronoiDiagram::compute_RVD()
	{
		face_begin = 0;
		face_end = p_Mesh->meshFacets.getFacetsNumber();

		Polygon F;
		/*
			compute the RVD of each facet and each seed
		*/
		for (t_index t = face_begin; t < face_end; ++t)
		{
			printf(" ------ facet %d --------\n", t);
			Facet temp_facet = p_Mesh->meshFacets.getFacet(t);
			current_center = computeCenter(temp_facet.m_v1, temp_facet.m_v2, temp_facet.m_v3);
			if (t == 185)
			{
				face_begin = face_begin; 
			}
			current_near_points = findNearestPoints(current_center, seeds_n);

			for (t_index i = 0; i < seeds_n; ++i)
			{
				current_seed = current_near_points[i];
				/*if (current_seed == 7)
				{
					current_seed++;
					current_seed--;
				}*/
				F.initialize_from_mesh_facet(p_Mesh, t);
				
				/*
					compute the intersection between a cell and a facet

					 current polygon : the clipped area with a cell to a facet
				*/
				current_polygon = intersect_cell_facet(current_seed, F);

				if (current_polygon->getVertex_nb() == 0)
					break;

				//muniplate the current polygon 
				polygonHandler->clear();
				polygonHandler->setPolygon(current_polygon);
				polygonHandler->compute_centriod();
				
				polygon_nb[current_seed]++;
				/*
					store the information of the polygon to the exact seed
				*/
				seedsUpdater.addInformation(polygonHandler->getCurrentPosition(), polygonHandler->getCurrentWeight(), current_seed);
			}
		}

		seedsUpdater.UpdateSeeds();

		p_Points->savePointsWithSeedStore(seedsUpdater.getSeedsPosition());

		return true;
	}

	Polygon* RestrictedVoronoiDiagram::intersect_cell_facet(t_index seed, Polygon& F)
	{
		seed_neighbors = findNearestPoints(seed, seeds_n);

		Polygon* ping = &F;
		Polygon* pong = &polygon_buffer;

		for (int i = 0; i < seeds_n; ++i)
		{
			int j = seed_neighbors[i];
			//test the clip_by_plane's correctness
			

			//printf("CPU :计算 %d 和 %d的半平面情况\n", seed, j);
			//printf("---------before---------\n");
			//printf("ping\n");
			//ping->show_polygon();
			//printf("pong\n");
			//pong->show_polygon();
			clip_by_plane(*ping, *pong, seed, (t_index)j);
			//swap ping and pong
			swap_polygons(ping, pong);
			//printf("---------after----------\n");
			//printf("ping\n");
			//ping->show_polygon();
			//printf("pong\n");
			//pong->show_polygon();
		}
		//swap the ping polygon and pong polygon
		//the result must be stored in the ping polygons
		return ping;
	}

	void RestrictedVoronoiDiagram::swap_polygons(Polygon*& ping, Polygon*& pong)
	{
		if (ping != &empty_polygon && ping != &polygon_buffer) {
			// First clipping operation, ping points to F
			// (current facet copied)
			ping = &polygon_buffer;
			pong = &empty_polygon;
		}
		else {
			P_RVD::geo_swap(ping, pong);
		}
	}

	void RestrictedVoronoiDiagram::clip_by_plane(Polygon& ping, Polygon& pong, t_index i, t_index j)
	{
		// pong is the result of ping
		ping.clip_by_plane(pong, *p_Points, i, j);
	}
}