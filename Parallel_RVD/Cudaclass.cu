#include "cuda_runtime.h"
#include "device_launch_parameters.h"

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
*/
__device__ double computeCentriod(double a, double b, double c, double centriod)
{
	
}

/*
	a thread to compute the RVD in parallel way
*/
__global__ void compute_RVD(double* seeds_pointer, int seeds_nb, 
							double* mesh_vertex,   int mesh_vertex_nb,
							int* mesh_facet,    int mesh_facet_nb)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	//
}


extern "C" void runCuda(double* host_seeds_pointer, double* host_mesh_vertex_pointer,
	int* host_facet_index, int points_nb, int mesh_vertex_nb, int mesh_facet_nb)
{
	double* dev_seeds_pointer;
	double* dev_mesh_vertex_pointer;
	int* dev_mesh_facet_index;

	checkCudaErrors(cudaSetDevice(0));

	//allocate the memory
	checkCudaErrors(cudaMalloc((void**)&dev_seeds_pointer, sizeof(double) * points_nb * 3));
	checkCudaErrors(cudaMalloc((void**)&dev_mesh_vertex_pointer, sizeof(double) * mesh_vertex_nb * 3));
	checkCudaErrors(cudaMalloc((void**)&dev_mesh_facet_index, sizeof(int) * mesh_facet_nb * 3));

	//pass the data from host to device
	checkCudaErrors(cudaMemcpy(dev_seeds_pointer, host_seeds_pointer, sizeof(double) * points_nb * 3, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_mesh_vertex_pointer, host_mesh_vertex_pointer, sizeof(double) * points_nb * 3, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_mesh_facet_index, host_facet_index, sizeof(int) * points_nb * 3, cudaMemcpyHostToDevice));

	int x = (int)sqrt(mesh_facet_nb);
	int y = (int)(mesh_facet_nb / x + 1);
	//set the block
	compute_RVD << < x, y >> >(dev_seeds_pointer, points_nb, dev_mesh_vertex_pointer, mesh_vertex_nb, dev_mesh_facet_index, mesh_facet_nb);

}