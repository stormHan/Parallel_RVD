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
void clip_by_plane(Cuda_Polygon& ping, Cuda_Polygon& pong, int seed_i, int seed_j, double* seeds_pointer, int seeds_nb)
{
	//reset the pong
	pong.vertex_nb = 0;

	if (ping.vertex_nb == 0)
		return;

	//get the position of seed i and seed j
	double3 position_i = { seeds_pointer[seed_i * 3 + 0], seeds_pointer[seed_i * 3 + 1], seeds_pointer[seed_i * 3 + 2] };
	double3 position_j = { seeds_pointer[seed_j * 3 + 0], seeds_pointer[seed_j * 3 + 1], seeds_pointer[seed_j * 3 + 2] };

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
				I.neigh_s = (seed_j);
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
__device__
void intersection_clip_facet_with_knn(Cuda_Polygon& current_polygon , int current_seed, double* seeds_pointer, int seeds_nb, int *seeds_neighbor_index, int k)
{
	
	//set a buffer pointer to store the polygon
	Cuda_Polygon polygon_buffer;
	
	for (int i = 0; i < k; ++i)
	{
		int j = seeds_neighbor_index[current_seed * k + i];
		if (current_seed != j)
		{
			clip_by_plane(current_polygon, polygon_buffer, current_seed, j, seeds_pointer, seeds_nb);
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
	if (triangle_nb > 0){

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

	double3 v1 = { mesh_vertex[facet_index.x * 3 + 0],
		mesh_vertex[facet_index.x * 3 + 1],
		mesh_vertex[facet_index.x * 3 + 2] };

	double3 v2 = { mesh_vertex[facet_index.y * 3 + 0],
		mesh_vertex[facet_index.y * 3 + 1],
		mesh_vertex[facet_index.y * 3 + 2] };

	double3 v3 = { mesh_vertex[facet_index.z * 3 + 0],
		mesh_vertex[facet_index.z * 3 + 1],
		mesh_vertex[facet_index.z * 3 + 2] };

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
	
	/*if (tid == 0)
	{
		ret_seeds[0] = current_polygon.vertex[0].x;
		ret_seeds[1] = current_polygon.vertex[0].y;
		ret_seeds[2] = current_polygon.vertex[0].z;
		ret_seeds[3] = current_polygon.vertex[0].w;
		ret_seeds[4] = current_polygon.vertex[0].neigh_s;
	}*/
	
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

		intersection_clip_facet_with_knn(current_polygon, current_seed, seeds_pointer, seeds_nb, seeds_neighbor_index, k);
		
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
	}

	//__syncthreads();
	
	//for (int i = 0; i < seeds_nb * 4; ++i)
	//{
	//	ret_seeds[i] = SeedsInformation[i];
	//}
	//for (int i = 0; i < seeds_nb; ++i)
	//{
	//	ret_seeds[i] = SeedsPolygon_nb[i];
	//}
	return;
}


/*
	a thread to compute the RVD in parallel way
*/
//__global__ 
//void compute_RVD(double* seeds_pointer, int seeds_nb,
//double* mesh_vertex, int mesh_vertex_nb,
//int* mesh_facet, int mesh_facet_nb, double* ret_seeds,
//double* test_dis, double* test_centriod,
//int* test_index, double* test_seeds,
//int* test_polygon_nb)
//{
//	int tid = blockIdx.x * blockDim.x + threadIdx.x;
//
//	//有多少个facet,开启多少个线程,用tid进行索引
//	if (tid < mesh_facet_nb){
//		//得到对应facet的三个索引
//		int f_idx1 = mesh_facet[tid * 3 + 0];
//		int f_idx2 = mesh_facet[tid * 3 + 1];
//		int f_idx3 = mesh_facet[tid * 3 + 2];
//
//		int3 facet_index = { f_idx1, f_idx2, f_idx3 };
//
//		double3 v1 = { mesh_vertex[facet_index.x * 3 + 0],
//			mesh_vertex[facet_index.x * 3 + 1],
//			mesh_vertex[facet_index.x * 3 + 2] };
//
//		double3 v2 = { mesh_vertex[facet_index.y * 3 + 0],
//			mesh_vertex[facet_index.y * 3 + 1],
//			mesh_vertex[facet_index.y * 3 + 2] };
//
//		double3 v3 = { mesh_vertex[facet_index.z * 3 + 0],
//			mesh_vertex[facet_index.z * 3 + 1],
//			mesh_vertex[facet_index.z * 3 + 2] };
//
//
//		double3 centriod = computeCentriod(v1, v2, v3);
//
//		
//		int* nearest_points = (int*)malloc(sizeof(int) * 20);
//
//		findNearestPoints(centriod, seeds_pointer, seeds_nb, nearest_points, 20);
//
//		if (tid == 0)
//		{
//			test_dis[0] = v1.x;
//			test_dis[1] = v1.y;
//			test_dis[2] = v1.z;
//			test_dis[3] = v2.x;
//			test_dis[4] = v2.y;
//			test_dis[5] = v2.z;
//			test_dis[6] = v3.x;
//			test_dis[7] = v3.y;
//			test_dis[8] = v3.z;
//
//			test_centriod[0] = centriod.x;
//			test_centriod[1] = centriod.y;
//			test_centriod[2] = centriod.z;
//
//			test_index[21] = seeds_nb;
//
//			for (int i = 0; i < 20; ++i)
//			{
//				test_index[i] = nearest_points[i];
//			}
//		}
//		return;
//
//		double* polygon = (double*)malloc(sizeof(double) * 10 * 4);
//		int vertex_nb;
//		/*
//			get the nearest 20 points' index
//			clip the facet
//			*/
//		for (int i = 0; i < 20; ++i)
//		{
//			int current_seed = nearest_points[i];
//
//			vertex_nb = 3;
//
//			//intialize the polygon with the 3 vertex of the current facet
//			/*
//				polygon pointer can be apart by several vertex
//				a vertex is made up with x,y,z,w
//				w : the weight of a vertex
//				*/
//			polygon[0] = v1.x; polygon[1] = v1.y; polygon[2] = v1.z; polygon[3] = 1.0;
//			polygon[4] = v2.x; polygon[5] = v2.y; polygon[6] = v2.z; polygon[7] = 1.0;
//			polygon[8] = v3.x; polygon[9] = v3.y; polygon[10] = v3.z; polygon[11] = 1.0;
//
//			intersection_clip_facet(polygon, vertex_nb, current_seed, seeds_pointer, seeds_nb);
//
//			//now we get the clipped polygon stored in "polygon"
//			//take care of the sychonize
//			//change the polygon data into "weight" and "position"
//			if (vertex_nb == 0)
//				break;
//
//			double weight;
//			double3 position;
//
//			int _v1 = 0;
//			int _v2, _v3;
//
//			double3 pos1, pos2, pos3;
//			double d1, d2, d3;
//			int triangle_nb = vertex_nb - 2;
//
//			double total_weight = 0.0;
//			double3 centriodTimesWeight = { 0.0, 0.0, 0.0 };
//
//			double current_weight = 0.0;
//			double3 current_posTimesWeight = { 0.0, 0.0, 0.0 };
//
//			for (int i = 1; i < vertex_nb - 1; ++i)
//			{
//				_v2 = i; _v3 = i + 1;
//
//				pos1 = { polygon[_v1 * 4 + 0], polygon[_v1 * 4 + 1], polygon[_v1 * 4 + 2] };
//				d1 = polygon[_v1 * 4 + 3];
//
//				pos2 = { polygon[_v2 * 4 + 0], polygon[_v2 * 4 + 1], polygon[_v2 * 4 + 2] };
//				d2 = polygon[_v2 * 4 + 3];
//
//				pos3 = { polygon[_v3 * 4 + 0], polygon[_v3 * 4 + 1], polygon[_v3 * 4 + 2] };
//				d3 = polygon[_v3 * 4 + 3];
//
//				computeTriangleCentriod(pos1, pos2, pos3, d1, d2, d3, centriodTimesWeight, total_weight);
//
//				current_weight += total_weight;
//				current_posTimesWeight.x += centriodTimesWeight.x;
//				current_posTimesWeight.y += centriodTimesWeight.y;
//				current_posTimesWeight.z += centriodTimesWeight.z;
//
//				total_weight = 0.0;
//				centriodTimesWeight = { 0.0, 0.0, 0.0 };
//			}
//
//			//atomicAdd(&SeedsPolygon_nb[current_seed], 1);
//			if (current_weight != 0 && triangle_nb > 0){
//				current_posTimesWeight.x /= current_weight;
//				current_posTimesWeight.y /= current_weight;
//				current_posTimesWeight.z /= current_weight;
//
//				current_weight /= triangle_nb;
//
//				double3 temp_pos;
//
//				temp_pos.x = current_posTimesWeight.x * current_weight;
//				temp_pos.y = current_posTimesWeight.y * current_weight;
//				temp_pos.z = current_posTimesWeight.z * current_weight;
//
//				MyAtomicAdd(&SeedsInformation[current_seed * 4 + 0], temp_pos.x);
//				MyAtomicAdd(&SeedsInformation[current_seed * 4 + 1], temp_pos.y);
//				MyAtomicAdd(&SeedsInformation[current_seed * 4 + 2], temp_pos.z);
//				MyAtomicAdd(&SeedsInformation[current_seed * 4 + 3], current_weight);
//
//
//			}
//			/*------------------------ test part -------------------------*/
//			if (tid == 0 && i == 0)
//			{
//				test_index[0] = current_seed;
//				for (int i = 0; i < 20; ++i)
//				{
//					test_index[i + 1] = nearest_points[i];
//				}
//
//				test_dis[0] = v1.x;
//				test_dis[1] = v1.y;
//				test_dis[2] = v1.z;
//				test_dis[3] = v2.x;
//				test_dis[4] = v2.y;
//				test_dis[5] = v2.z;
//				test_dis[6] = v3.x;
//				test_dis[7] = v3.y;
//				test_dis[8] = v3.z;
//			}
//
//		}
//		__syncthreads();
//		for (int i = 0; i < 100; ++i)
//		{
//			ret_seeds[i] = SeedsInformation[i];
//		}
//		free(nearest_points);
//		free(polygon);
//	}
//}
//
//extern "C" void runCuda(double* host_seeds_pointer, double* host_mesh_vertex_pointer,
//	int* host_facet_index, int points_nb, int mesh_vertex_nb, int mesh_facet_nb)
//{
//	double* dev_seeds_pointer;
//	double* dev_mesh_vertex_pointer;
//	int* dev_mesh_facet_index;
//
//	double* dev_seedsInformation;
//
//	double* host_test_seeds = (double*)malloc(sizeof(double) *  points_nb * 4);
//	double* dev_test_seeds;
//	//test data
//	
//	double* host_test_dis = (double*)malloc(sizeof(double) * points_nb);
//	double* dev_test_dis;
//
//	double* host_test_cen = (double*)malloc(sizeof(double) * 3);
//	double* dev_test_center;
//
//	int* dev_seedsPolyonNumber;
//
//	int* host_test_index = (int*)malloc(sizeof(int) * points_nb);
//	int* dev_test_index;
//
//	int* host_test_polygon = (int*)malloc(sizeof(int) * points_nb);
//	int* dev_test_polygon;
//	
//	checkCudaErrors(cudaSetDevice(0));
//	
//	//allocate the memory
//	checkCudaErrors(cudaMalloc((void**)&dev_seeds_pointer, sizeof(double) * points_nb * 3));
//	checkCudaErrors(cudaMalloc((void**)&dev_mesh_vertex_pointer, sizeof(double) * mesh_vertex_nb * 3));
// 	checkCudaErrors(cudaMalloc((void**)&dev_mesh_facet_index, sizeof(int) * mesh_facet_nb * 3));
// 	checkCudaErrors(cudaMalloc((void**)&dev_seedsInformation, sizeof(double) * points_nb * 4));
//	//checkCudaErrors(cudaMalloc((void**)&dev_seedsPolyonNumber, sizeof(int) * points_nb));
//	checkCudaErrors(cudaMalloc((void**)&dev_test_seeds, sizeof(double) * points_nb * 4));
//
//
//	//allocate test memory
//
//	checkCudaErrors(cudaMalloc((void**)&dev_test_dis, sizeof(double) * points_nb));
//	checkCudaErrors(cudaMalloc((void**)&dev_test_center, sizeof(double) * 3));
//	checkCudaErrors(cudaMalloc((void**)&dev_test_index, sizeof(int) * points_nb));
//	
//	checkCudaErrors(cudaMalloc((void**)&dev_test_polygon, sizeof(int) * points_nb));
//	
//	//pass the data from host to device
//	checkCudaErrors(cudaMemcpy(dev_seeds_pointer, host_seeds_pointer, sizeof(double) * points_nb * 3, cudaMemcpyHostToDevice));
//	checkCudaErrors(cudaMemcpy(dev_mesh_vertex_pointer, host_mesh_vertex_pointer, sizeof(double) * mesh_vertex_nb * 3, cudaMemcpyHostToDevice));
//	checkCudaErrors(cudaMemcpy(dev_mesh_facet_index, host_facet_index, sizeof(int) * mesh_facet_nb * 3, cudaMemcpyHostToDevice));
//
//	checkCudaErrors(cudaMemcpyToSymbol(SeedsInformation, &dev_seedsInformation, sizeof(double*), size_t(0), cudaMemcpyHostToDevice));
//	/*checkCudaErrors(cudaMemcpyToSymbol(SeedsPolygon_nb, &dev_seedsPolyonNumber, sizeof(int*), size_t(0), cudaMemcpyHostToDevice));*/
//	
//	int x = (int)sqrt(mesh_facet_nb);
//	int y = (int)(mesh_facet_nb / x + 1);
//	//set the block
//	compute_RVD << < x, y >> >(dev_seeds_pointer, points_nb, dev_mesh_vertex_pointer, mesh_vertex_nb, dev_mesh_facet_index, mesh_facet_nb, dev_test_seeds,
//		dev_test_dis, dev_test_center, dev_test_index, dev_test_seeds, dev_test_polygon);
//	
//	checkCudaErrors(cudaDeviceSynchronize());
//
//	checkCudaErrors(cudaMemcpy(host_test_dis, dev_test_dis, sizeof(double) * points_nb, cudaMemcpyDeviceToHost));
//	checkCudaErrors(cudaMemcpy(host_test_cen, dev_test_center, sizeof(double) * 3, cudaMemcpyDeviceToHost));
//	checkCudaErrors(cudaMemcpy(host_test_index, dev_test_index, sizeof(int) * points_nb, cudaMemcpyDeviceToHost));
//	checkCudaErrors(cudaMemcpy(host_test_polygon, dev_test_polygon, sizeof(int) * points_nb, cudaMemcpyDeviceToHost));
//	
//	checkCudaErrors(cudaMemcpy(host_test_seeds, dev_test_seeds, sizeof(double) * points_nb * 4, cudaMemcpyDeviceToHost));
//	printf("Thread 1:\n");
//	printf("centriod : %.6lf, %.6lf, %.6lf \n", host_test_cen[0], host_test_cen[1], host_test_cen[2]);
//	printf("seeds_nb : %d \n", host_test_index[21]);
//	printf("v1 : %.6lf, %.6lf, %.6lf \n", host_test_dis[0], host_test_dis[1], host_test_dis[2]);
//	printf("v2 : %.6lf, %.6lf, %.6lf \n", host_test_dis[3], host_test_dis[4], host_test_dis[5]);
//	printf("v3 : %.6lf, %.6lf, %.6lf \n", host_test_dis[6], host_test_dis[7], host_test_dis[8]);
//	printf("nearest 20 points' index:\n");
//	for (int i = 0; i < 20; ++i)
//	{
//		printf("%d  : %d  \n",i, host_test_index[i]);
//	}
//	return;
//	//printf("current seed : %d\n", host_test_index[41]);
//	//printf("current seed position : %.17lf, %.17lf, %.17lf \n", host_test_dis[3], host_test_dis[4], host_test_dis[5]);
//
//	//printf("the 20 seeds near the current seed\n");
//	//for (int i = 0; i < 20; ++i)
//	//{
//	//	printf("%d  : %d  \n",i, host_test_index[i + 21]);
//	//}
//
//	//printf(" %d 和 %d 确定平面\n", host_test_index[41], host_test_index[42]);
//	//printf("原polygon:\n");
//	//printf("vertex number : %d\n", host_test_index[43]);
//
//	//for (int i = 0; i < host_test_index[43]; ++i)
//	//{
//	//	printf("x : %.17lf, y: %.17lf, z : %.17lf, w : %.6lf\n\n", host_test_dis[6 + i * 4 + 0], host_test_dis[6 + i * 4 + 1], host_test_dis[6 + i * 4 + 2], host_test_dis[6 + i * 4 + 3]);
//	//}
//	//printf("现polygon:\n");
//	//printf("vertex number : %d   buffer vertex number : %d\n", host_test_index[44], host_test_index[45]);
//	//for (int i = 0; i < host_test_index[44]; ++i)
//	//{
//	//	printf("x : %.17lf, y: %.17lf, z : %.17lf, w : %.6lf\n\n", host_test_dis[101 + i * 4 + 0], host_test_dis[101 + i * 4 + 1], host_test_dis[101 + i * 4 + 2], host_test_dis[101 + i * 4 + 3]);
//	//}
//
//	//printf("增加的点\n");
//	//for (int i = 0; i < 4; ++i)
//	//{
//	//	printf("x : %.17lf, y: %.17lf, z : %.17lf, w : %.6lf\n\n", host_test_dis[150 + i * 4 + 0], host_test_dis[150 + i * 4 + 1], host_test_dis[150 + i * 4 + 2], host_test_dis[150 + i * 4 + 3]);
//	//}
//
//	//printf("thread 1, facet 1, current seed : %d clipped polygon vertex number : %d\n",host_test_index[20], host_test_index[88]);
//	//printf("successfully running!\n");
//
//	ofstream fout("d.txt");
//	//printf("test : %d current seed :%d  vertex number : %d\n",host_test_index[1], host_test_index[2], host_test_index[3]);
//	//fout << "seed 185  nearest points" << endl;
//	//for (int i = 0; i < 20; ++i)
//	//{
//	//	fout << host_test_index[i] << endl;
//	//}
//	////fout << "v1 : " << host_test_dis[0] << host_test_dis[1] << host_test_dis[2] << endl;
//	////fout << "v2 : " << host_test_dis[0] << host_test_dis[1] << host_test_dis[2] << endl;
//	////fout << "v2 : " << host_test_dis[0] << host_test_dis[1] << host_test_dis[2] << endl;
//
//	for (int i = 0; i < points_nb; ++i)
//	{
//		if (host_test_seeds[i * 4 + 3] != 0){
//			host_test_seeds[i * 4 + 0] /= host_test_seeds[i * 4 + 3];
//			host_test_seeds[i * 4 + 1] /= host_test_seeds[i * 4 + 3];
//			host_test_seeds[i * 4 + 2] /= host_test_seeds[i * 4 + 3];
//		}
//		//printf("seeds %d : %.16lf, %.16lf, %.16lf, , %.16lf \n", i, host_test_seeds[i * 4 + 0]
//			//, host_test_seeds[i * 4 + 1], host_test_seeds[i * 4 + 2], host_test_seeds[i * 4 + 3]);
//		
//		fout << i << ':' << setprecision(16) <<host_test_seeds[i * 4 + 0] << ' ' << host_test_seeds[i * 4 + 1] << ' ' << host_test_seeds[i * 4 + 2] << ' ' << host_test_seeds[i * 4 + 3] << endl;
//	}
//	/*for (int i = 0; i < points_nb; ++i)
//	{
//		fout << i << ':' << host_test_polygon[i] << endl;
//	}*/
//}

extern "C" void runRVD(double* host_seeds_pointer, double* host_mesh_vertex_pointer,
	int* host_facet_index, int points_nb, int mesh_vertex_nb, int mesh_facet_nb,
	std::vector<int> facet_center_neigbors, std::vector<int> seeds_neighbors, std::vector<int>& seeds_polygon_nb)
{
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

	int threads = 512;
	int blocks = mesh_facet_nb / threads + ((mesh_facet_nb % threads) ? 1 : 0);
	
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	float elapsed_time;
	cudaEventRecord(start, 0);
	compute_RVD_with_knn << <threads, blocks >> >(dev_seeds_pointer, points_nb, dev_mesh_vertex_pointer, mesh_vertex_nb, dev_mesh_facet_index, mesh_facet_nb, dev_facet_center_neighbors,
		dev_seeds_neighbors, 10, ret_seeds, test_seeds);

	CheckCUDAError("kenerl function");
	
	cudaMemcpy(host_seedsInfo, ret_seeds, sizeof(double) * points_nb * 4, cudaMemcpyDeviceToHost);

	CheckCUDAError("pass data back to CPU");
	
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsed_time, start, stop);
	printf("Compute RVD time: %lfms\n", elapsed_time);
	getchar();
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
	cudaFree(dev_facet_center_neighbors);
	cudaFree(dev_seeds_neighbors);

	return;
}