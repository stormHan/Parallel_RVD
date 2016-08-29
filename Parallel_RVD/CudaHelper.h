/*
	header to help use the CUDA
	some basics works
*/

#ifndef H_CUDAHELPER_H
#define H_CUDAHELPER_H

#include "Points.h"



namespace P_RVD
{
	/*
	help trans the points data

	input:
	p : the input points data
	output:
	host_data_pointer: pointer in cpu
	device_data_pointer: pointer in gpu
	*/
	void trans_points(const Points p, double* host_data_pointer)
	{
		int number = p.getPointsNumber();
		host_data_pointer = (double*)malloc(3 * sizeof(double) * number);

		for (t_index t = 0; t < number; ++t)
		{
			host_data_pointer[3 * t + 0] = p.getPoint(t).x;
			host_data_pointer[3 * t + 1] = p.getPoint(t).y;
			host_data_pointer[3 * t + 2] = p.getPoint(t).z;
		}
	}

	/*
		help trans the mesh data

		input:
		m : the input mesh
		output:
		host_data_pointer: the mesh vertices data(double)
		host_index_pointer: the mesh facet indice(int)
	*/
	void trans_mesh(const Mesh m, double* host_data_pointer, int* host_index_pointer)
	{
		int vertex_nb = m.meshVertices.getPointNumber();
		host_data_pointer = (double*)malloc(3 * sizeof(double) * vertex_nb);

		for (t_index t = 0; t < vertex_nb; ++t)
		{
			host_data_pointer[3 * t + 0] = m.meshVertices.getPoint(t).x;
			host_data_pointer[3 * t + 1] = m.meshVertices.getPoint(t).y;
			host_data_pointer[3 * t + 2] = m.meshVertices.getPoint(t).z;
		}

		int facet_nb = m.meshFacets.getFacetsNumber();
		host_index_pointer = (int*)malloc(3 * sizeof(int) * facet_nb);

		for (t_index t = 0; t < facet_nb; ++t)
		{
			host_index_pointer[3 * t + 0] = m.meshFacets.getFacet(t).m_v1;
			host_index_pointer[3 * t + 1] = m.meshFacets.getFacet(t).m_v2;
			host_index_pointer[3 * t + 2] = m.meshFacets.getFacet(t).m_v3;
		}
	}
}


#endif /* H_CUDAHELPER_H */