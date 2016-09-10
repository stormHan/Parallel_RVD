#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "Math_basics.cuh"

#include "Polygon.h"
#include <iostream>

using namespace P_RVD;
/*
	process the cuda error
*/
inline void checkCudaErrors(cudaError err)
{
	if (cudaSuccess != err)
	{
		fprintf(stderr, "CUDA Runtime API error : %s.\n", cudaGetErrorString(err));
		return;
	}
}

/*
	device function
	conmpute 3 vertex's centroid of a facet
*/
__device__ double3 computeCentriod(double3 a, double3 b, double3 c)
{
	double3 ret = { (a.x + b.x + c.x) / 3,
		(a.y + b.y + c.y) / 3,
		(a.z + b.z + c.z) / 3, };
	return ret;
}

/*
	device function
	find nearest n seeds in device

	input:
	centriod : the physical position of the seed
	seeds_pointer : pointer which store the position of the seeds
	seeds_nb : the number of seeds

	output:
	nearest_points : the result of the function which store the index of the nearest seeds
*/
__device__
void findNearestPoints(double3 centriod, double* seeds_pointer, int seeds_nb, int* nearest_points, int nearest_points_nb)
{
	double* distance_from_centroid = (double*)malloc(sizeof(double) * seeds_nb);

	//find the nearest points
	for (int i = 0; i < seeds_nb; ++i)
	{
		double3 seeds_pos = { seeds_pointer[i * 3 + 0], seeds_pointer[i * 3 + 1], seeds_pointer[i * 3 + 2] };
		distance_from_centroid[i] = computeDistance(centriod, seeds_pos);
	}

	//设置一个数组来存储下标
	int* index = (int*)malloc(sizeof(int) * seeds_nb);
	for (int i = 0; i < seeds_nb; ++i)
	{
		index[i] = i;
	}

	for (int i = 1; i < seeds_nb; ++i)
	{
		double key = distance_from_centroid[i];
		int store_index = index[i];
		for (int j = i - 1; j >= 0; --j)
		{
			if (distance_from_centroid[j] >= key)
			{
				distance_from_centroid[j + 1] = distance_from_centroid[j];
				index[j + 1] = index[j];
				if (j == 0)
				{
					distance_from_centroid[j] = key;
					index[j] = store_index;
				}
			}
			else
			{
				distance_from_centroid[j + 1] = key;
				index[j + 1] = store_index;
				break;
			}

		}
	}
	for (int i = 0; i < 20; ++i)
	{
		nearest_points[i] = index[i];
	}
}


/*
	device function
	compute the intersection of a facet and a seed(a cell)

	input:

*/
__device__
void intersection_clip_facet(double* polygon_pointer, int vertex_nb, int current_seed, double* seeds_pointer, int seeds_nb)
{
	double3 current_seed_pos = { seeds_pointer[current_seed * 3 + 0],
		seeds_pointer[current_seed * 3 + 1],
		seeds_pointer[current_seed * 3 + 2] };

	int* seeds_neighbors = (int*)malloc(sizeof(int) * 20);
	findNearestPoints(current_seed_pos, seeds_pointer, seeds_nb, seeds_neighbors, 20);

	//set a buffer pointer to store the polygon
	double* polygon_buffer = (double*)malloc(sizeof(double) * 10 * 4);
	int polygon_buffer_nb = 0;

	for (int i = 0; i < 20; ++i)
	{
		int j = seeds_neighbors[i];

		clip_by_plane(polygon_pointer, vertex_nb, polygon_buffer, polygon_buffer_nb, current_seed, j, seeds_pointer, seeds_nb);
		swap_polygons(polygon_pointer, vertex_nb, polygon_buffer, polygon_buffer_nb);
	}
}

/*
	device function
	
	compute the intersection(stored in a polygon) clipped by a bisector defined by seed i and j

	input : polygon ping and its number
	output: polygon pong and its number
*/
__device__
void clip_by_plane(double* ping, int ping_nb, double* pong, int pong_nb, int seed_i, int seed_j, double* seeds_pointer, int seeds_nb)
{
	//reset the pong
	pong_nb = 0;

	if (ping_nb == 0)
		return;

	//get the position of seed i and seed j
	double3 position_i = { seeds_pointer[seed_i * 3 + 0], seeds_pointer[seed_i * 3 + 1], seeds_pointer[seed_i * 3 + 2] };
	double3 position_j = { seeds_pointer[seed_j * 3 + 0], seeds_pointer[seed_j * 3 + 1], seeds_pointer[seed_j * 3 + 2] };

	// Compute d = n . (2m), where n is the
	// normal vector of the bisector [i, j]
	// and m the middle point of the bisector.
	double d;
	d = dot(position_i,position_j)

}

/*
	device function
	swap the ping and pong to make sure 
	ping store the result;
*/
__device__
void swap_polygons(double* ping, int ping_nb, double* pong, int pong_nb)
{
	//!!! can be accerate
	//不确定是否可以使用memset memcpy
	double temp;
	int temp_nb = ping_nb;

	ping_nb = pong_nb;
	pong_nb = temp_nb;

	for (int i = 0; i < 40; ++i)
	{
		temp = ping[i];
		ping[i] = pong[i];
		pong[i] = temp;
	}
}


/*
	a thread to compute the RVD in parallel way
*/
__global__ void compute_RVD(double* seeds_pointer, int seeds_nb, 
							double* mesh_vertex,   int mesh_vertex_nb,
							int* mesh_facet,    int mesh_facet_nb,
							double* test_dis, double* test_centriod,
							int* test_index)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	//有多少个facet,开启多少个线程,用tid进行索引
	if (tid < mesh_facet_nb){
		//得到对应facet的三个索引
		int f_idx1 = mesh_facet[tid * 3 + 0];
		int f_idx2 = mesh_facet[tid * 3 + 1];
		int f_idx3 = mesh_facet[tid * 3 + 2];

		int3 facet_index = { f_idx1, f_idx2, f_idx3 };

		double3 v1 = { mesh_vertex[facet_index.x * 3 + 0],
			mesh_vertex[facet_index.x * 3 + 1],
			mesh_vertex[facet_index.x * 3 + 2] };

		double3 v2 = { mesh_vertex[facet_index.y * 3 + 0],
			mesh_vertex[facet_index.y * 3 + 1],
			mesh_vertex[facet_index.y * 3 + 2] };

		double3 v3 = { mesh_vertex[facet_index.z * 3 + 0],
			mesh_vertex[facet_index.z * 3 + 1],
			mesh_vertex[facet_index.z * 3 + 2] };


		double3 centriod = computeCentriod(v1, v2, v3);

		int* nearest_points = (int*)malloc(sizeof(int) * 20);

		findNearestPoints(centriod, seeds_pointer, seeds_nb, nearest_points, 20);
		
		/*
			get the nearest 20 points' index
			clip the facet
		*/
		for (int i = 0; i < 20; ++i)
		{
			int current_seed = nearest_points[i];

			double* polygon = (double*)malloc(sizeof(double) * 10 * 4);
			int vertex_nb = 3;
			//intialize the polygon with the 3 vertex of the current facet
			/*
				polygon pointer can be apart by several vertex
				a vertex is made up with x,y,z,w
				w : the weight of a vertex
			*/
			polygon[0] = v1.x; polygon[1] = v1.y; polygon[2] = v1.z; polygon[3] = 1.0;
			polygon[4] = v2.x; polygon[5] = v2.y; polygon[6] = v2.z; polygon[7] = 1.0;
			polygon[8] = v3.x; polygon[9] = v3.y; polygon[10] = v3.z; polygon[11] = 1.0;


		}
		

		//pass the data to test data
		if (tid == 0)
		{
			test_centriod[0] = centriod.x;
			test_centriod[1] = centriod.y;
			test_centriod[2] = centriod.z;

			for (int i = 0; i < 20; ++i)
			{
				test_index[i] = nearest_points[i];
			}
		}
	}
}


extern "C" void runCuda(double* host_seeds_pointer, double* host_mesh_vertex_pointer,
	int* host_facet_index, int points_nb, int mesh_vertex_nb, int mesh_facet_nb)
{
	double* dev_seeds_pointer;
	double* dev_mesh_vertex_pointer;
	int* dev_mesh_facet_index;

	//test data
	double* host_test_dis = (double*)malloc(sizeof(double) * points_nb);
	double* dev_test_dis;

	double* host_test_cen = (double*)malloc(sizeof(double) * 3);
	double* dev_test_center;

	int* host_test_index = (int*)malloc(sizeof(int) * points_nb);
	int* dev_test_index;


	checkCudaErrors(cudaSetDevice(0));
	
	//allocate the memory
	checkCudaErrors(cudaMalloc((void**)&dev_seeds_pointer, sizeof(double) * points_nb * 3));
	checkCudaErrors(cudaMalloc((void**)&dev_mesh_vertex_pointer, sizeof(double) * mesh_vertex_nb * 3));
	checkCudaErrors(cudaMalloc((void**)&dev_mesh_facet_index, sizeof(int) * mesh_facet_nb * 3));


	//allocate test memory
	checkCudaErrors(cudaMalloc((void**)&dev_test_dis, sizeof(double) * points_nb));
	checkCudaErrors(cudaMalloc((void**)&dev_test_center, sizeof(double) * 3));
	checkCudaErrors(cudaMalloc((void**)&dev_test_index, sizeof(int) * points_nb));

	//pass the data from host to device
	checkCudaErrors(cudaMemcpy(dev_seeds_pointer, host_seeds_pointer, sizeof(double) * points_nb * 3, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_mesh_vertex_pointer, host_mesh_vertex_pointer, sizeof(double) * points_nb * 3, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_mesh_facet_index, host_facet_index, sizeof(int) * points_nb * 3, cudaMemcpyHostToDevice));

	int x = (int)sqrt(mesh_facet_nb);
	int y = (int)(mesh_facet_nb / x + 1);
	//set the block
	compute_RVD << < x, y >> >(dev_seeds_pointer, points_nb, dev_mesh_vertex_pointer, mesh_vertex_nb, dev_mesh_facet_index, mesh_facet_nb,
		dev_test_dis, dev_test_center, dev_test_index);

	checkCudaErrors(cudaMemcpy(host_test_dis, dev_test_dis, sizeof(double) * points_nb, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(host_test_cen, dev_test_center, sizeof(double) * 3, cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(host_test_index, dev_test_index, sizeof(int) * points_nb, cudaMemcpyDeviceToHost));

	printf("center position : %.17lf, %.17lf, %.17lf \n", host_test_cen[0], host_test_cen[1], host_test_cen[2]);
	for (int i = 0; i < points_nb; ++i)
	{
		printf("%d  : %.16lf  \n",i, host_test_dis[i]);
	}

	for (int i = 0; i < 20; ++i)
	{
		printf("%d  : %d  \n",i, host_test_index[i]);
	}

	printf("successfully running!");
}