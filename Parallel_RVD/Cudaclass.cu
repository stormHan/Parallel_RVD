#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "Math_basics.cuh"

#include <iostream>
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

		//pass the data to test data
		if (tid == 0)
		{
			test_centriod[0] = centriod.x;
			test_centriod[1] = centriod.y;
			test_centriod[2] = centriod.z;

			for (int i = 0; i < seeds_nb; ++i)
			{
				test_dis[i] = distance_from_centroid[i];
				test_index[i] = index[i];
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

	for (int i = 0; i < points_nb; ++i)
	{
		printf("%d  : %d  \n",i, host_test_index[i]);
	}

	printf("successfully running!");
}