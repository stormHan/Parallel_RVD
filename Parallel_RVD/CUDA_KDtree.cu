#include "CUDA_KDtree.h"
#include <cuda_runtime.h>
#include <cstdio>

#define CUDA_STACK 1000

namespace P_RVD{

	void CheckCUDAError(const char *msg)
	{
		cudaError_t err = cudaGetLastError();
		if (cudaSuccess != err) {
			fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString(err));
			exit(EXIT_FAILURE);
		}
	}
	__device__ double Distance(const Kd_tree_point &a, const Kd_tree_point &b)
	{
		double dist = 0;

		for (int i = 0; i < KDTREE_DIM; i++) {
			double d = a.coords[i] - b.coords[i];
			dist += d*d;
		}

		return dist;
	}

	__device__ void SearchAtNode(const CUDA_KDNode *nodes, const int *indexes, const Kd_tree_point *pts, int cur, const Kd_tree_point &query, int *ret_index, double *ret_dist, int *ret_node)
	{
		// Finds the first potential candidate

		int best_idx = 0;
		double best_dist = FLT_MAX;

		while (true) {
			int split_axis = nodes[cur].level % KDTREE_DIM;

			if (nodes[cur].left == -1) {
				*ret_node = cur;

				for (int i = 0; i < nodes[cur].num_indexes; i++) {
					int idx = indexes[nodes[cur].indexes + i];
					double dist = Distance(query, pts[idx]);
					if (dist < best_dist) {
						best_dist = dist;
						best_idx = idx;
					}
				}

				break;
			}
			else if (query.coords[split_axis] < nodes[cur].split_value) {
				cur = nodes[cur].left;
			}
			else {
				cur = nodes[cur].right;
			}
		}

		*ret_index = best_idx;
		*ret_dist = best_dist;
	}


	__device__ void cudaSearchAtNodeRange(const CUDA_KDNode *nodes, const int *indexes, const Kd_tree_point *pts, const Kd_tree_point &query, int cur, double range, int *ret_index, double *ret_dist)
	{
		// Goes through all the nodes that are within "range"

		int best_idx = 0;
		double best_dist = FLT_MAX;

		// Ok, we don't have nice STL vectors to use, and we can't dynamically allocate memory with CUDA??
		// We'll use a fixed length stack, increase this as required
		int to_visit[CUDA_STACK];
		int to_visit_pos = 0;

		to_visit[to_visit_pos++] = cur;

		while (to_visit_pos) {
			int next_search[CUDA_STACK];
			int next_search_pos = 0;

			while (to_visit_pos) {
				cur = to_visit[to_visit_pos - 1];
				to_visit_pos--;

				int split_axis = nodes[cur].level % KDTREE_DIM;

				if (nodes[cur].left == -1) {
					for (int i = 0; i < nodes[cur].num_indexes; i++) {
						int idx = indexes[nodes[cur].indexes + i];
						double d = Distance(query, pts[idx]);

						if (d < best_dist) {
							best_dist = d;
							best_idx = idx;
						}
					}
				}
				else {
					double d = query.coords[split_axis] - nodes[cur].split_value;

					// There are 3 possible scenarios
					// The hypercircle only intersects the left region
					// The hypercircle only intersects the right region
					// The hypercricle intersects both

					if (fabs(d) > range) {
						if (d < 0)
							next_search[next_search_pos++] = nodes[cur].left;
						else
							next_search[next_search_pos++] = nodes[cur].right;
					}
					else {
						next_search[next_search_pos++] = nodes[cur].left;
						next_search[next_search_pos++] = nodes[cur].right;
					}
				}
			}

			// No memcpy available??
			for (int i = 0; i < next_search_pos; i++)
				to_visit[i] = next_search[i];

			to_visit_pos = next_search_pos;
		}

		*ret_index = best_idx;
		*ret_dist = best_dist;
	}


	__device__ void Search(const CUDA_KDNode *nodes, const int *indexes, const Kd_tree_point *pts, const Kd_tree_point &query, int *ret_index, double *ret_dist)
	{
		// Find the first closest node, this will be the upper bound for the next searches
		int best_node = 0;
		int best_idx = 0;
		double best_dist = FLT_MAX;
		double radius = 0;

		SearchAtNode(nodes, indexes, pts, 0 /* root */, query, &best_idx, &best_dist, &best_node);

		radius = sqrt(best_dist);

		// Now find other possible candidates
		int cur = best_node;

		while (nodes[cur].parent != -1) {
			// Go up
			int parent = nodes[cur].parent;
			int split_axis = nodes[parent].level % KDTREE_DIM;

			// Search the other node
			double tmp_dist = FLT_MAX;
			int tmp_idx;

			if (fabs(nodes[parent].split_value - query.coords[split_axis]) <= radius) {
				// Search opposite node
				if (nodes[parent].left != cur)
					cudaSearchAtNodeRange(nodes, indexes, pts, query, nodes[parent].left, radius, &tmp_idx, &tmp_dist);
				else
					cudaSearchAtNodeRange(nodes, indexes, pts, query, nodes[parent].right, radius, &tmp_idx, &tmp_dist);
			}

			if (tmp_dist < best_dist) {
				best_dist = tmp_dist;
				best_idx = tmp_idx;
			}

			cur = parent;
		}

		*ret_index = best_idx;
		*ret_dist = best_dist;
	}

	__device__ void SearchAtNode_knn(const CUDA_KDNode *nodes, const int *indexes, const Kd_tree_point *pts, int cur, const Kd_tree_point &query, int *ret_index, double *ret_dist, int *ret_node, int k)
	{
		int neighbor_nb = 0;


		//travasel to the nodes
		while (true)
		{
			int split_axis = nodes[cur].level % KDTREE_DIM;

			if (nodes[cur].left == -1){
				//Get to the leaf node
				
				neighbor_nb += nodes[cur].num_indexes;
				while (neighbor_nb < k)
				{
					cur = nodes[cur].parent;
					neighbor_nb = nodes[cur].num_indexes;
				}
				*ret_node = cur;
				
				//Now we get enough neighbors in cur node
				//int *temp_index = (int*)malloc(sizeof(int) * nodes[cur].num_indexes);
				//double *temp_dists = (double*)malloc(sizeof(double) * nodes[cur].num_indexes);
				int temp_index[20];
				double temp_dists[20];
				
				for (int i = 0; i < nodes[cur].num_indexes; ++i)
				{
					temp_index[i] = indexes[nodes[cur].indexes + i];
					temp_dists[i] = Distance(query, pts[temp_index[i]]);
				}
				
				int n = nodes[cur].num_indexes;
				//利用k次冒泡得到前小的距离
				int best_idx = 0;
				double best_dist = FLT_MAX;
				for (int i = 0; i < k; ++i)
				{
					for (int j = i; j < n; ++j)
					{
						if (temp_dists[j] < best_dist)
						{
							best_dist = temp_dists[j];
							best_idx = temp_index[j];

							temp_dists[j] = temp_dists[i];
							temp_index[j] = temp_index[i];

							temp_index[i] = best_idx;
							temp_dists[i] = best_dist;
						}
					}

					ret_index[i] = best_idx;
					ret_dist[i] = best_dist;

					best_idx = 0;
					best_dist = FLT_MAX;
				}
				//free(temp_index);
				//free(temp_dists);
				break;
			}
			else if (query.coords[split_axis] < nodes[cur].split_value){
				cur = nodes[cur].left;
			}
			else{
				cur = nodes[cur].right;
			}
		}

	}

	__device__ void SearchAtiNodeRange_knn(const CUDA_KDNode *nodes, const int *indexes, const Kd_tree_point *pts, const Kd_tree_point &query, int cur, double range, int* &ret_index, double* &ret_dist, int k)
	{
		// Goes through all the nodes that within "radius"
		int to_visit[CUDA_STACK];
		int to_visit_pos = 0;

		to_visit[to_visit_pos++] = cur;

		while (to_visit_pos){
			int next_search[CUDA_STACK];
			int next_search_pos = 0;

			while (to_visit_pos){
				cur = to_visit[to_visit_pos - 1];
				to_visit_pos--;

				int split_axis = nodes[cur].level % KDTREE_DIM;

				//already get to the leaf
				if (nodes[cur].left == -1){
					for (int i = 0; i < nodes[cur].num_indexes; ++i)
					{
						int idx = indexes[nodes[cur].indexes + i];
						double d = Distance(query, pts[idx]);

						for (int j = 0; j < k; ++j)
						{

							if (d < ret_dist[j])
							{
								for (int q = k - 1; q > j; --q)
								{
									ret_index[q] = ret_index[q - 1];
									ret_dist[q] = ret_dist[q - 1];
								}
								ret_index[j] = idx;
								ret_dist[j] = d;
								break;
							}
						}
					}
				}
				else{
					double d = query.coords[split_axis] - nodes[cur].split_value;

					// There are 3 possible scenarios
					// The hypercircle only intersects the left region
					// The hypercircle only intersects the right region
					// The hypercricle intersects both

					if (fabs(d) > range) {
						if (d < 0)
							next_search[next_search_pos++] = nodes[cur].left;
						else
							next_search[next_search_pos++] = nodes[cur].right;
					}
					else {
						next_search[next_search_pos++] = nodes[cur].left;
						next_search[next_search_pos++] = nodes[cur].right;
					}
				}
			}
			// No memcpy available??
			for (int i = 0; i < next_search_pos; i++)
				to_visit[i] = next_search[i];

			to_visit_pos = next_search_pos;
		}
	}


	__device__ void Search_knn(const CUDA_KDNode *nodes, const int *indexes, const Kd_tree_point *pts, const Kd_tree_point &query, int *ret_index, double *ret_dist, int k)
	{
		// Find the first closest node, this will be the upper bound for the next searches
		int best_node = 0;
		int* k_idx = (int*)malloc(sizeof(int) * k);
		double* k_dist = (double*)malloc(sizeof(double) * k);
		double radius = 0;

		SearchAtNode_knn(nodes, indexes, pts, 0 /* root */, query, k_idx, k_dist, &best_node, k);

		radius = sqrt(k_dist[k - 1]);

		//Now find other posiible candidates
		int cur = best_node;

		while (nodes[cur].parent != -1)
		{
			int parent = nodes[cur].parent;
			int split_axis = nodes[parent].level % KDTREE_DIM;

			//Search the other node

			if (fabs(nodes[parent].split_value - query.coords[split_axis]) <= radius)
			{
				
				if (nodes[parent].left != cur)
					SearchAtiNodeRange_knn(nodes, indexes, pts, query, nodes[parent].left, radius, k_idx, k_dist, k);
				else
					SearchAtiNodeRange_knn(nodes, indexes, pts, query, nodes[parent].right, radius, k_idx, k_dist, k);
					
			}
			cur = parent;
		}

		for (int i = 0; i < k; ++i)
		{
			ret_index[i] = k_idx[i];
			ret_dist[i] = k_dist[i];
		}
		free(k_idx);
		free(k_dist);
	}


	__global__ void SearchBatch_knn(const CUDA_KDNode *nodes, const int *indexes, const Kd_tree_point *pts, int num_pts, Kd_tree_point *queries, int num_queries, int *ret_index, double *ret_dist, int k)
	{
		int idx = blockIdx.x*blockDim.x + threadIdx.x;

		if (idx >= num_queries)
			return;

		//test(nodes, indexes, pts, queries[idx], &ret_index[idx * k], &ret_dist[idx * k], k);
		Search_knn(nodes, indexes, pts, queries[idx], &ret_index[idx * k], &ret_dist[idx * k], k);
	}

	__global__ void SearchBatch(const CUDA_KDNode *nodes, const int *indexes, const Kd_tree_point *pts, int num_pts, Kd_tree_point *queries, int num_queries, int *ret_index, double *ret_dist)
	{
		int idx = blockIdx.x*blockDim.x + threadIdx.x;

		if (idx >= num_queries)
			return;

		Search(nodes, indexes, pts, queries[idx], &ret_index[idx], &ret_dist[idx]);
	}

	CUDA_KDTree::~CUDA_KDTree()
	{
		cudaFree(m_gpu_indexes);
		cudaFree(m_gpu_nodes);
		cudaFree(m_gpu_points);
	}

	void CUDA_KDTree::CreateKDtree(KDNode *root, int num_nodes, const std::vector<Kd_tree_point> &data, int max_level)
	{
		m_num_points = data.size();

		cudaMalloc((void**)&m_gpu_nodes, sizeof(CUDA_KDNode) * num_nodes);
		cudaMalloc((void**)&m_gpu_indexes, sizeof(int) * m_num_points * (max_level + 1));
		cudaMalloc((void**)&m_gpu_points, sizeof(Kd_tree_point) * m_num_points);

		CheckCUDAError("CreateKDtree");

		std::vector <CUDA_KDNode> cpu_nodes(num_nodes);
		std::vector <int> indexes(m_num_points * (max_level + 1));
		std::vector <KDNode*> to_visit;

		int cur_pos = 0;
		to_visit.push_back(root);

		while (to_visit.size())
		{
			std::vector<KDNode*> next_search;

			while (to_visit.size())
			{
				KDNode *cur = to_visit.back();
				to_visit.pop_back();

				int id = cur->id;
				cpu_nodes[id].level = cur->level;
				cpu_nodes[id].parent = cur->_parent;
				cpu_nodes[id].left = cur->_left;
				cpu_nodes[id].right = cur->_right;
				cpu_nodes[id].split_value = cur->split_value;
				cpu_nodes[id].num_indexes = cur->index.size();

				if (cur->index.size()) {
					for (unsigned int i = 0; i < cur->index.size(); i++)
						indexes[cur_pos + i] = cur->index[i];

					cpu_nodes[id].indexes = cur_pos;
					cur_pos += cur->index.size();
				}
				else {
					cpu_nodes[id].indexes = -1;
				}

				if (cur->left)
					next_search.push_back(cur->left);

				if (cur->right)
					next_search.push_back(cur->right);
			}

			to_visit = next_search;
		}
		cudaMemcpy(m_gpu_nodes, &cpu_nodes[0], sizeof(CUDA_KDNode) * cpu_nodes.size(), cudaMemcpyHostToDevice);
		cudaMemcpy(m_gpu_indexes, &indexes[0], sizeof(int) * indexes.size(), cudaMemcpyHostToDevice);
		cudaMemcpy(m_gpu_points, &data[0], sizeof(Kd_tree_point) * data.size(), cudaMemcpyHostToDevice);

		CheckCUDAError("CreateKDTree");
	}

	void CUDA_KDTree::Search(const std::vector<Kd_tree_point> &queries, std::vector<int> &indexes, std::vector<double> &dists)
	{
		int threads = 512;
		int blocks = queries.size() / threads + ((queries.size() % threads) ? 1 : 0);

		Kd_tree_point *gpu_queries;
		int *gpu_ret_indexes;
		double *gpu_ret_dist;

		indexes.resize(queries.size());
		dists.resize(queries.size());

		cudaMalloc((void**)&gpu_queries, sizeof(double) * queries.size() * KDTREE_DIM);
		cudaMalloc((void**)&gpu_ret_indexes, sizeof(int) * queries.size() * KDTREE_DIM);
		cudaMalloc((void**)&gpu_ret_dist, sizeof(double) * queries.size() * KDTREE_DIM);

		CheckCUDAError("Search pass 1");

		cudaMemcpy(gpu_queries, &queries[0], sizeof(double) * queries.size() * KDTREE_DIM, cudaMemcpyHostToDevice);

		CheckCUDAError("Search pass 2");

		printf("CUDA blocks/threads: %d %d\n", blocks, threads);

		SearchBatch << <blocks, threads >> >(m_gpu_nodes, m_gpu_indexes, m_gpu_points, m_num_points, gpu_queries, queries.size(), gpu_ret_indexes, gpu_ret_dist);
		cudaThreadSynchronize();

		CheckCUDAError("Search");

		cudaMemcpy(&indexes[0], gpu_ret_indexes, sizeof(int)*queries.size(), cudaMemcpyDeviceToHost);
		cudaMemcpy(&dists[0], gpu_ret_dist, sizeof(double)*queries.size(), cudaMemcpyDeviceToHost);

		cudaFree(gpu_queries);
		cudaFree(gpu_ret_indexes);
		cudaFree(gpu_ret_dist);
	}

	void CUDA_KDTree::Search(const Mesh &query_mesh, std::vector<int> &indexes, std::vector<double> &dists)
	{
		int threads = 512;
		int blocks = query_mesh.meshFacets.getFacetsNumber() / threads + ((query_mesh.meshFacets.getFacetsNumber() % threads) ? 1 : 0);

		std::vector<Kd_tree_point> queries(query_mesh.meshFacets.getFacetsNumber());
		
		//algorithm 1 compute the centriod of a triangle in cpu
		for (int i = 0; i < query_mesh.meshFacets.getFacetsNumber(); ++i)
		{
			Facet temp_f = query_mesh.meshFacets.getFacet(i);
			
			Vector3d p1 = query_mesh.meshVertices.getPoint(temp_f.m_v1);
			Vector3d p2 = query_mesh.meshVertices.getPoint(temp_f.m_v2);
			Vector3d p3 = query_mesh.meshVertices.getPoint(temp_f.m_v3);

			Vector3d center = Math::computeCenter(p1, p2, p3);
			queries[i].coords[0] = center.x;
			queries[i].coords[1] = center.y;
			queries[i].coords[2] = center.z;
		}

		Kd_tree_point *gpu_queries;
		int *gpu_ret_indexes;
		double *gpu_ret_dist;

		indexes.resize(queries.size());
		dists.resize(queries.size());

		cudaMalloc((void**)&gpu_queries, sizeof(double) * queries.size() * KDTREE_DIM);
		cudaMalloc((void**)&gpu_ret_indexes, sizeof(int) * queries.size() * KDTREE_DIM);
		cudaMalloc((void**)&gpu_ret_dist, sizeof(double) * queries.size() * KDTREE_DIM);

		CheckCUDAError("Search pass 1");

		cudaMemcpy(gpu_queries, &queries[0], sizeof(double) * queries.size() * KDTREE_DIM, cudaMemcpyHostToDevice);

		CheckCUDAError("Search pass 2");

		printf("CUDA blocks/threads: %d %d\n", blocks, threads);

		SearchBatch << <blocks, threads >> >(m_gpu_nodes, m_gpu_indexes, m_gpu_points, m_num_points, gpu_queries, queries.size(), gpu_ret_indexes, gpu_ret_dist);
		cudaThreadSynchronize();

		CheckCUDAError("Search");

		cudaMemcpy(&indexes[0], gpu_ret_indexes, sizeof(int)*queries.size(), cudaMemcpyDeviceToHost);
		cudaMemcpy(&dists[0], gpu_ret_dist, sizeof(double)*queries.size(), cudaMemcpyDeviceToHost);

		cudaFree(gpu_queries);
		cudaFree(gpu_ret_indexes);
		cudaFree(gpu_ret_dist);

	}

	void CUDA_KDTree::Search_knn(const std::vector<Kd_tree_point> &queries, std::vector<int> &indexes, std::vector<double> &dists, int k)
	{
		int threads = 512;
		int blocks = queries.size() / threads + ((queries.size() % threads) ? 1 : 0);

		Kd_tree_point *gpu_queries;
		int *gpu_ret_indexes;
		double *gpu_ret_dist;

		indexes.resize(queries.size() * k);
		dists.resize(queries.size() * k);

		cudaMalloc((void**)&gpu_queries, sizeof(double) * queries.size() * KDTREE_DIM);
		cudaMalloc((void**)&gpu_ret_indexes, sizeof(int) * k * queries.size());
		cudaMalloc((void**)&gpu_ret_dist, sizeof(double) * k * queries.size());
		CheckCUDAError("Initialize the gpu pointer");

		//copy the query data
		cudaMemcpy(gpu_queries, &queries[0], sizeof(double) * queries.size() * KDTREE_DIM, cudaMemcpyHostToDevice);
		CheckCUDAError("Copy the data");

		// Variables for duration evaluation
		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		float elapsed_time;

		cudaEventRecord(start, 0);
			SearchBatch_knn << < blocks, threads >> >(m_gpu_nodes, m_gpu_indexes, m_gpu_points, m_num_points, gpu_queries, queries.size(), gpu_ret_indexes, gpu_ret_dist, k);
		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&elapsed_time, start, stop);
		printf("Searching %d points, My kdtree done in %f ms \n",queries.size(), elapsed_time);
		// Destroy cuda event
		cudaEventDestroy(start);
		cudaEventDestroy(stop);

		
		//cudaThreadSynchronize();

		CheckCUDAError("kernel function");

		//Copy back the data from GPU
		cudaMemcpy(&indexes[0], gpu_ret_indexes, sizeof(int) * queries.size() * k, cudaMemcpyDeviceToHost);
		cudaMemcpy(&dists[0], gpu_ret_dist, sizeof(double) * k * queries.size(), cudaMemcpyDeviceToHost);
		CheckCUDAError("Copy back the data from GPU");

		cudaFree(gpu_queries);
		cudaFree(gpu_ret_dist);
		cudaFree(gpu_ret_indexes);

	}

	void CUDA_KDTree::Search_knn(const Mesh &query_mesh, std::vector<int> &indexes, std::vector<double> &dists, int k)
	{
		int threads = 512;
		int blocks = query_mesh.meshFacets.getFacetsNumber() / threads + ((query_mesh.meshFacets.getFacetsNumber() % threads) ? 1 : 0);

		std::vector<Kd_tree_point> queries(query_mesh.meshFacets.getFacetsNumber());
		//algorithm 1 compute the centriod of a triangle in cpu
		for (int i = 0; i < query_mesh.meshFacets.getFacetsNumber(); ++i)
		{
			Facet temp_f = query_mesh.meshFacets.getFacet(i);

			Vector3d p1 = query_mesh.meshVertices.getPoint(temp_f.m_v1);
			Vector3d p2 = query_mesh.meshVertices.getPoint(temp_f.m_v2);
			Vector3d p3 = query_mesh.meshVertices.getPoint(temp_f.m_v3);

			Vector3d center = Math::computeCenter(p1, p2, p3);
			queries[i].coords[0] = center.x;
			queries[i].coords[1] = center.y;
			queries[i].coords[2] = center.z;
		}

		Kd_tree_point *gpu_queries;
		int *gpu_ret_indexes;
		double *gpu_ret_dist;

		indexes.resize(queries.size() * k);
		dists.resize(queries.size() * k);

		cudaMalloc((void**)&gpu_queries, sizeof(double) * queries.size() * KDTREE_DIM);
		cudaMalloc((void**)&gpu_ret_indexes, sizeof(int) * k * queries.size());
		cudaMalloc((void**)&gpu_ret_dist, sizeof(double) * k * queries.size());
		CheckCUDAError("Initialize the gpu pointer");

		//copy the query data
		cudaMemcpy(gpu_queries, &queries[0], sizeof(double) * queries.size() * KDTREE_DIM, cudaMemcpyHostToDevice);
		CheckCUDAError("Copy the data");

		// Variables for duration evaluation
		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		float elapsed_time;

		cudaEventRecord(start, 0);
		SearchBatch_knn << < blocks, threads >> >(m_gpu_nodes, m_gpu_indexes, m_gpu_points, m_num_points, gpu_queries, queries.size(), gpu_ret_indexes, gpu_ret_dist, k);
		cudaEventRecord(stop, 0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&elapsed_time, start, stop);
		printf("Searching %d points, My kdtree done in %f ms \n", queries.size(), elapsed_time);
		// Destroy cuda event
		cudaEventDestroy(start);
		cudaEventDestroy(stop);

		
		//cudaThreadSynchronize();

		CheckCUDAError("kernel function");

		//Copy back the data from GPU
		cudaMemcpy(&indexes[0], gpu_ret_indexes, sizeof(int) * queries.size() * k, cudaMemcpyDeviceToHost);
		cudaMemcpy(&dists[0], gpu_ret_dist, sizeof(double) * k * queries.size(), cudaMemcpyDeviceToHost);
		CheckCUDAError("Copy back the data from GPU");

		cudaFree(gpu_queries);
		cudaFree(gpu_ret_dist);
		cudaFree(gpu_ret_indexes);
	}
}