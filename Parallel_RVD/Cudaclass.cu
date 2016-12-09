#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "device_atomic_functions.h"

#include "Math_basics.cuh"
#include "AtomicAction.cuh"

#include "Polygon.h"
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace P_RVD;

//set global variable to store the seeds' position and weight
__device__ double* SeedsInformation;

__device__ int* SeedsPolygon_nb;

const int CUDA_Stack_size = 10;

texture<int2, 1> t_vertex;
texture<int2, 1> t_points;
texture<int, 1> t_points_nn;

struct Cuda_Vertex
{
	double x;
	double y;
	double z;
	double w;

	int neigh_s = -1;
};

struct Cuda_Polygon
{
	Cuda_Vertex vertex[10];
	int vertex_nb;
};
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

void CheckCUDAError(const char *msg)
{
	cudaError_t err = cudaGetLastError();
	if (cudaSuccess != err) {
		fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}

}
__device__
double MyAtomicAdd(double* address, double val)
{
	unsigned long long int* address_as_ull = (unsigned long long int*)address;

	unsigned long long int old = *address_as_ull, assumed;

	do{
		assumed = old;
		old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
	} while (assumed != old);

	return __longlong_as_double(old);
}

static __inline__ __device__
double fetch_double(texture<int2, 1> t, int i){
	int2 v = tex1Dfetch(t, i);
	return __hiloint2double(v.y, v.x);
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

compute the intersection(stored in a polygon) clipped by a bisector defined by seed i and j

input : polygon ping and its number
output: polygon pong and its number
*/
__device__
void clip_by_plane(Cuda_Polygon& ping, Cuda_Polygon& pong, double3 position_i,
double3 position_j, int j)
{
	//reset the pong
	pong.vertex_nb = 0;

	if (ping.vertex_nb == 0)
		return;

	// Compute d = n . (2m), where n is the
	// normal vector of the bisector [i, j]
	// and m the middle point of the bisector.
	double d = 0.0;
	d = dot(add(position_i, position_j), sub(position_i, position_j));

	//The predecessor of the first vertex is the last vertex
	int prev_k = ping.vertex_nb - 1;
	
	//get the position data
	Cuda_Vertex* prev_vk = &ping.vertex[prev_k];

	double3 prev_vertex_position = { prev_vk->x, prev_vk->y, prev_vk->z };

	//then we compute prev_vertex_position "cross" n 
	//prev_l = prev_vertex_position . n

	double prev_l = dot(prev_vertex_position, sub(position_i, position_j));

	int prev_status = sgn(2.0 * prev_l - d);

	//traverse the Vertex in this Polygon
	for (int k = 0; k < ping.vertex_nb; ++k)
	{
		Cuda_Vertex* vk = &ping.vertex[k];
		double3 vertex_position = { vk->x, vk->y, vk->z };

		double l = dot(vertex_position, sub(position_i, position_j));
		int status = sgn(2.0 * l - d);

		//If status of edge extremities differ,
		//then there is an intersection.
		if (status != prev_status && (prev_status) != 0)
		{
			// create the intersection and update the Polyon
			Cuda_Vertex I;

			//compute the position and weight
			double denom = 2.0 * (prev_l - l);
			double lambda1, lambda2;

			// Shit happens!
			if (m_fabs(denom) < 1e-20)
			{
				lambda1 = 0.5;
				lambda2 = 0.5;
			}
			else
			{
				lambda1 = (d - 2.0 * l) / denom;
				// Note: lambda2 is also given
				// by (2.0*l2-d)/denom
				// (but 1.0 - lambda1 is a bit
				//  faster to compute...)
				lambda2 = 1.0 - lambda1;
			}

			//Set the Position of Vertex
			I.x = lambda1 * prev_vertex_position.x + lambda2 * vertex_position.x;
			I.y = lambda1 * prev_vertex_position.y + lambda2 * vertex_position.y;
			I.z = lambda1 * prev_vertex_position.z + lambda2 * vertex_position.z;

			//Set the Weight of Vertex
			I.w = (lambda1 * prev_vk->w + lambda2 * vk->w);

			if (status > 0)
			{
				I.neigh_s = (j);
			}
			else {
				I.neigh_s = (vk->neigh_s);
			}

			//add I to pong
			pong.vertex[pong.vertex_nb] = I;
			pong.vertex_nb++;

		}

		if (status > 0)
		{
			//add vertex to pong
			pong.vertex[pong.vertex_nb] = *vk;
			pong.vertex_nb++;
		}

		prev_vk = vk;
		prev_vertex_position = vertex_position;
		prev_status = status;
		prev_l = l;
		prev_k = k;
	}
	return;
}

/*
device function
swap the ping and pong to make sure
ping store the result;
*/
__device__
void swap_polygons(Cuda_Polygon& ping, Cuda_Polygon& pong)
{
	//!!! can be accerate
	//不确定是否可以使用memset memcpy
	Cuda_Polygon t = ping;
	ping = pong;
	pong = t;
}

/*
*/
//__device__
//void intersection_clip_facet_with_knn(Cuda_Polygon& current_polygon , int current_seed, double* seeds_pointer, int seeds_nb, int *seeds_neighbor_index, int k)
//{
//	
//	//set a buffer pointer to store the polygon
//	Cuda_Polygon polygon_buffer;
//	
//	for (int i = 0; i < k; ++i)
//	{
//		int j = seeds_neighbor_index[current_seed * k + i];
//		if (current_seed != j)
//		{
//			clip_by_plane(current_polygon, polygon_buffer, current_seed, j, seeds_pointer, seeds_nb);
//			swap_polygons(current_polygon, polygon_buffer);
//		}
//	}
//	return;
//}

/*
*/
__device__
void intersection_clip_facet_SR(Cuda_Polygon& current_polygon, int i, double* seeds_pointer, int seeds_nb,
int *seeds_neighbor_index, int k)
{
	Cuda_Polygon polygon_buffer;

	double3 pi = {
		fetch_double(t_points, i * 3 + 0),
		fetch_double(t_points, i * 3 + 1),
		fetch_double(t_points, i * 3 + 2)
	};

	for (int t = 0; t < k; ++t)
	{
		int j = tex1Dfetch(t_points_nn, i * k + t);
		if (i != j)
		{
			double3 pj = {
				fetch_double(t_points, j * 3 + 0),
				fetch_double(t_points, j * 3 + 1),
				fetch_double(t_points, j * 3 + 2)
			};

			double dij =  distance2(pi, pj);
			double R2 = 0.0;
			
			for (int ii = 0; ii < current_polygon.vertex_nb; ++ii)
			{
				double3 pk = { current_polygon.vertex[ii].x, current_polygon.vertex[ii].y, current_polygon.vertex[ii].z };
				double dik = distance2(pi, pk);
				R2 = max(R2, dik);
			}
			if (dij > 4.1 * R2)
			{
				return;
			}

			clip_by_plane(current_polygon, polygon_buffer, pi, pj, j);
			swap_polygons(current_polygon, polygon_buffer);
		}
	}
	return;
}

__device__
void action(const Cuda_Polygon polygon, int current_seed)
{
	double weight;
	double3 position;

	int _v1 = 0;
	int _v2, _v3;

	double3 pos1, pos2, pos3;
	double d1, d2, d3;
	int triangle_nb = polygon.vertex_nb - 2;

	double total_weight = 0.0;
	double3 centriodTimesWeight = { 0.0, 0.0, 0.0 };

	double current_weight = 0.0;
	double3 current_posTimesWeight = { 0.0, 0.0, 0.0 };

	for (int i = 1; i < polygon.vertex_nb - 1; ++i)
	{
		_v2 = i; _v3 = i + 1;

		pos1 = { polygon.vertex[_v1].x, polygon.vertex[_v1].y, polygon.vertex[_v1].z };
		d1 = polygon.vertex[_v1].w;

		pos2 = { polygon.vertex[_v2].x, polygon.vertex[_v2].y, polygon.vertex[_v2].z };
		d2 = polygon.vertex[_v2].w;

		pos3 = { polygon.vertex[_v3].x, polygon.vertex[_v3].y, polygon.vertex[_v3].z };
		d3 = polygon.vertex[_v3].w;

		computeTriangleCentriod(pos1, pos2, pos3, d1, d2, d3, centriodTimesWeight, total_weight);

		current_weight += total_weight;
		current_posTimesWeight.x += centriodTimesWeight.x;
		current_posTimesWeight.y += centriodTimesWeight.y;
		current_posTimesWeight.z += centriodTimesWeight.z;

		total_weight = 0.0;
		centriodTimesWeight = { 0.0, 0.0, 0.0 };
	}
	atomicAdd(&SeedsPolygon_nb[current_seed], 1);
	//atomicAdd(&SeedsPolygon_nb[current_seed], 1);
	if (triangle_nb > 0){
		//atomicAdd(&SeedsPolygon_nb[current_seed], 1);
		current_weight /= triangle_nb;

		double3 temp_pos;

		temp_pos.x = current_posTimesWeight.x / triangle_nb;
		temp_pos.y = current_posTimesWeight.y / triangle_nb;
		temp_pos.z = current_posTimesWeight.z / triangle_nb;


		//try not to use the MyAtomicAdd
		//SeedsInformation[current_seed * 4 + 0] += temp_pos.x;
		//SeedsInformation[current_seed * 4 + 1] += temp_pos.y;
		//SeedsInformation[current_seed * 4 + 2] += temp_pos.z;
		//SeedsInformation[current_seed * 4 + 3] += current_weight;

		MyAtomicAdd(&SeedsInformation[current_seed * 4 + 0], temp_pos.x);
		MyAtomicAdd(&SeedsInformation[current_seed * 4 + 1], temp_pos.y);
		MyAtomicAdd(&SeedsInformation[current_seed * 4 + 2], temp_pos.z);
		MyAtomicAdd(&SeedsInformation[current_seed * 4 + 3], current_weight);
	}
}

__global__
void compute_RVD_with_knn(double* seeds_pointer, int seeds_nb,
double* mesh_vertex, int mesh_vertex_nb,
int* mesh_facet, int mesh_facet_nb,
int* facet_center_neighbor_index, int* seeds_neighbor_index, int k, double* ret_seeds, double* test_seeds)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= mesh_facet_nb) return;

	if (tid >= 0 && tid < seeds_nb){
		SeedsInformation[tid * 4 + 0] = 0.0;
		SeedsInformation[tid * 4 + 1] = 0.0;
		SeedsInformation[tid * 4 + 2] = 0.0;
		SeedsInformation[tid * 4 + 3] = 0.0;
	}

	int f_idx1 = mesh_facet[tid * 3 + 0];
	int f_idx2 = mesh_facet[tid * 3 + 1];
	int f_idx3 = mesh_facet[tid * 3 + 2];

	int3 facet_index = { f_idx1, f_idx2, f_idx3 };

	double3 v1 = {
		fetch_double(t_vertex, facet_index.x * 3 + 0),
		fetch_double(t_vertex, facet_index.x * 3 + 1),
		fetch_double(t_vertex, facet_index.x * 3 + 2)
	};
	double3 v2 = {
		fetch_double(t_vertex, facet_index.y * 3 + 0),
		fetch_double(t_vertex, facet_index.y * 3 + 1),
		fetch_double(t_vertex, facet_index.y * 3 + 2)
	};
	double3 v3 = {
		fetch_double(t_vertex, facet_index.z * 3 + 0),
		fetch_double(t_vertex, facet_index.z * 3 + 1),
		fetch_double(t_vertex, facet_index.z * 3 + 2)
	};

	Cuda_Polygon current_polygon;
	current_polygon.vertex_nb = 3;

	//intialize the polygon with the 3 vertex of the current facet
	/*
	polygon pointer can be apart by several vertex
	a vertex is made up with x,y,z,w
	w : the weight of a vertex
	*/
	current_polygon.vertex[0].x = v1.x; current_polygon.vertex[0].y = v1.y; current_polygon.vertex[0].z = v1.z; current_polygon.vertex[0].w = 1.0;
	current_polygon.vertex[1].x = v2.x; current_polygon.vertex[1].y = v2.y; current_polygon.vertex[1].z = v2.z; current_polygon.vertex[1].w = 1.0;
	current_polygon.vertex[2].x = v3.x; current_polygon.vertex[2].y = v3.y; current_polygon.vertex[2].z = v3.z; current_polygon.vertex[2].w = 1.0;
	Cuda_Polygon store = current_polygon;
	//doesn't have the stack?
	int to_visit[CUDA_Stack_size];
	int to_visit_pos = 0;

	int has_visited[CUDA_Stack_size];
	int has_visited_nb = 0;
	bool has_visited_flag = false;
	
	to_visit[to_visit_pos++] = facet_center_neighbor_index[tid];
	has_visited[has_visited_nb++] = to_visit[0];
	
	while (to_visit_pos){
		int current_seed = to_visit[to_visit_pos - 1];
		to_visit_pos--;

		intersection_clip_facet_SR(current_polygon, current_seed, seeds_pointer, seeds_nb, seeds_neighbor_index, k);
		/*
		if (tid == 0){
			ret_seeds[0] = current_polygon.vertex_nb;
			for (int i = 0; i < ret_seeds[0]; ++i){
				ret_seeds[1 + 5 * i + 0] = current_polygon.vertex[i].x;
				ret_seeds[1 + 5 * i + 1] = current_polygon.vertex[i].y;
				ret_seeds[1 + 5 * i + 2] = current_polygon.vertex[i].z;
				ret_seeds[1 + 5 * i + 3] = current_polygon.vertex[i].w;
				ret_seeds[1 + 5 * i + 4] = current_polygon.vertex[i].neigh_s;
			}
		}
		return;*/
		//use the RVD
		//action(current_polygon, current_seed);

		//now we get the clipped polygon stored in "polygon"
		//take care of the sychonize
		//change the polygon data into "weight" and "position"
		
		//Propagate to adjacent seeds
		for (int v = 0; v < current_polygon.vertex_nb; ++v)
		{
			Cuda_Vertex ve = current_polygon.vertex[v];
			int ns = ve.neigh_s;
			if (ns != -1)
			{
				for (int ii = 0; ii < has_visited_nb; ++ii)
				{
					//if the neighbor seed has clipped the polygon
					//the flag should be set "true"
					if (has_visited[ii] == ns)
						has_visited_flag = true;
				}
				//the neighbor seed is new!
				if (!has_visited_flag)
				{
					to_visit[to_visit_pos++] = ns;
					has_visited[has_visited_nb++] = ns;
				}
				has_visited_flag = false;
			}
		}
		current_polygon = store;
	}

	//__syncthreads();
	
	
	

	/*for (int i = 0; i < seeds_nb * 4; ++i)
	{
		ret_seeds[i] = seedsinformation[i];
	}*/
	/*for (int i = 0; i < seeds_nb; ++i)
	{
		ret_seeds[i] = SeedsPolygon_nb[i];
	}*/
	return;
}


extern "C" void runRVD(double* host_seeds_pointer, double* host_mesh_vertex_pointer,
	int* host_facet_index, int points_nb, int mesh_vertex_nb, int mesh_facet_nb,
	std::vector<int> facet_center_neigbors, std::vector<int> seeds_neighbors, std::vector<int>& seeds_polygon_nb)
{
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float elapsed_time;
	cudaEventRecord(start, 0);

	//GPU data
	double* dev_seeds_pointer;
	double* dev_mesh_vertex_pointer;
	int* dev_mesh_facet_index;
	double* dev_seedsInformation;
	double* ret_seeds;
	
	seeds_polygon_nb.resize(points_nb);
	int *dev_facet_center_neighbors, *dev_seeds_neighbors;
	
	/*---------test----------*/
	double* test_seeds;
	int* dev_seedsPolygonNumber;
	//CPU data
	double* host_seedsInfo = (double*)malloc(sizeof(double) * points_nb * 16);

	//allocate the memory
	checkCudaErrors(cudaMalloc((void**)&dev_seeds_pointer, sizeof(double) * points_nb * 3));
	checkCudaErrors(cudaMalloc((void**)&dev_mesh_vertex_pointer, sizeof(double) * mesh_vertex_nb * 3));
	checkCudaErrors(cudaMalloc((void**)&dev_mesh_facet_index, sizeof(int) * mesh_facet_nb * 3));
	checkCudaErrors(cudaMalloc((void**)&dev_seedsInformation, sizeof(double) * points_nb * 4));
	checkCudaErrors(cudaMalloc((void**)&ret_seeds, sizeof(double) * points_nb * 4));
	checkCudaErrors(cudaMalloc((void**)&test_seeds, sizeof(double) * points_nb * 16));

	checkCudaErrors(cudaMalloc((void**)&dev_seedsPolygonNumber, sizeof(int) * points_nb));

	checkCudaErrors(cudaMalloc((void**)&dev_facet_center_neighbors, sizeof(int) * facet_center_neigbors.size()));
	checkCudaErrors(cudaMalloc((void**)&dev_seeds_neighbors, sizeof(int) * seeds_neighbors.size()));

	checkCudaErrors(cudaMemcpyToSymbol(SeedsInformation, &dev_seedsInformation, sizeof(double*), size_t(0), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpyToSymbol(SeedsPolygon_nb, &dev_seedsPolygonNumber, sizeof(int*), size_t(0), cudaMemcpyHostToDevice));

	//pass the data from CPU to GPU
	cudaMemcpy(dev_seeds_pointer, host_seeds_pointer, sizeof(double) * points_nb * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_mesh_vertex_pointer, host_mesh_vertex_pointer, sizeof(double) * mesh_vertex_nb * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_mesh_facet_index, host_facet_index, sizeof(int) * mesh_facet_nb * 3, cudaMemcpyHostToDevice);
	cudaMemcpy(dev_facet_center_neighbors, &facet_center_neigbors[0], sizeof(int) * facet_center_neigbors.size(), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_seeds_neighbors, &seeds_neighbors[0], sizeof(int) * seeds_neighbors.size(), cudaMemcpyHostToDevice);

	CheckCUDAError("cudaMemcpyHostToDevice");

	cudaBindTexture(NULL, t_vertex, dev_mesh_vertex_pointer, sizeof(double) * mesh_vertex_nb * 3);
	cudaBindTexture(NULL, t_points, dev_seeds_pointer, sizeof(double) * points_nb * 3);
	cudaBindTexture(NULL, t_points_nn, dev_seeds_neighbors, sizeof(int) * seeds_neighbors.size());

	int threads = 512;
	int blocks = mesh_facet_nb / threads + ((mesh_facet_nb % threads) ? 1 : 0);
	
	
	compute_RVD_with_knn << <threads, blocks >> >(dev_seeds_pointer, points_nb, dev_mesh_vertex_pointer, mesh_vertex_nb, dev_mesh_facet_index, mesh_facet_nb, dev_facet_center_neighbors,
		dev_seeds_neighbors, 20, ret_seeds, test_seeds);

	CheckCUDAError("kenerl function");
	
	cudaMemcpy(host_seedsInfo, ret_seeds, sizeof(double) * points_nb * 4, cudaMemcpyDeviceToHost);

	CheckCUDAError("pass data back to CPU");
	
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventSynchronize(start);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsed_time, start, stop);
	printf("Compute RVD time: %lfms\n", elapsed_time);
	//getchar();
	//for (int i = 0; i < points_nb; ++i)
	//{
	//	if (fabs(host_seedsInfo[i * 4 + 3]) >= 1e-12){
	//		host_seedsInfo[i * 4 + 0] /= host_seedsInfo[i * 4 + 3];
	//		host_seedsInfo[i * 4 + 1] /= host_seedsInfo[i * 4 + 3];
	//		host_seedsInfo[i * 4 + 2] /= host_seedsInfo[i * 4 + 3];
	//	}
	//}
	//std::ofstream fileout("S2new.txt");
	//for (int i = 0; i < points_nb ; ++i)
	//{
	//	//printf("Line %d : x : %.17lf, y : %.17lf, z : %.17lf, w : %.17lf\n", i, host_seedsInfo[i * 4 + 0], host_seedsInfo[i * 4 + 1], host_seedsInfo[i * 4 + 2], host_seedsInfo[i * 4 + 3]);
	//	fileout << "Line " << i << ':' << setprecision(16) << host_seedsInfo[i * 4 + 0] << ' ' << host_seedsInfo[i * 4 + 1] << ' ' << host_seedsInfo[i * 4 + 2] << ' ' << host_seedsInfo[i * 4 + 3] << endl;
	//	//fileout << "Points " << i << ':' << host_seedsInfo[i] << endl;
	//	//seeds_polygon_nb[i] = host_seedsInfo[i];
	//}
	
	free(host_seedsInfo);
	cudaFree(dev_seeds_pointer);
	cudaFree(dev_mesh_vertex_pointer);
	cudaFree(dev_mesh_facet_index);
	cudaFree(dev_seedsInformation);
	cudaFree(ret_seeds);
	cudaFree(test_seeds);
	cudaFree(dev_facet_center_neighbors);
	cudaFree(dev_seeds_neighbors);

	cudaUnbindTexture(t_vertex);
	cudaUnbindTexture(t_points);
	cudaUnbindTexture(t_points_nn);

	return;
}