#include "CUDA_KDtree.h"
#include <cuda_runtime.h>
#include <cstdio>


namespace P_RVD{

	void CheckCUDAError(const char *msg)
	{
		cudaError_t err = cudaGetLastError();
		if (cudaSuccess != err) {
			fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString(err));
			exit(EXIT_FAILURE);
		}
	}

	CUDA_KDTree::~CUDA_KDTree()
	{
		cudaFree(m_gpu_indexes);
		cudaFree(m_gpu_nodes);
		cudaFree(m_gpu_points);
	}

	void CUDA_KDTree::CreateKDtree(KDNode *root, int num_nodes, const std::vector<Kd_tree_point> &data)
	{
		m_num_points = data.size();

		cudaMalloc((void**)&m_gpu_nodes, sizeof(CUDA_KDNode) * num_nodes);
		cudaMalloc((void**)&m_gpu_indexes, sizeof(int) * m_num_points);
		cudaMalloc((void**)&m_gpu_points, sizeof(Kd_tree_point) * m_num_points);

		CheckCUDAError("CreateKDtree");

		std::vector <CUDA_KDNode> cpu_nodes(num_nodes);
		std::vector <int> indexes(m_num_points);
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

		
	}
}