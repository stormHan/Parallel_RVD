/*
	CUDA kdtree
	used after we have created a cpu kdtree
*/

#ifndef H_CUDA_KDTREE_H
#define H_CUDA_KDTREE_H

#include "Kdtree.h"
#include "Mesh.h"
#include "Math_basics.h"
#include <vector>

namespace P_RVD
{
	/*
		cuda kdtree node

		every node has an unique number(index)

		level : indicates the height of this node, we can get the axis of split value by a node's level.
		parent, left, right : the number(index) of the node's parent, left child and right child.
		split value : the exact value to differ the points of the left child from the right.
		num_indexes
	*/
	struct CUDA_KDNode
	{
		int level;
		int parent, left, right;
		double split_value;
		int num_indexes;
		int indexes;
	};

	class CUDA_KDTree
	{
	public:
		~CUDA_KDTree();
		void CreateKDtree(KDNode *root, int num_nodes, const std::vector<Kd_tree_point> &data, int max_level);
		
		/*
			the search data : queries
		*/
		void Search(const std::vector<Kd_tree_point> &queries, std::vector<int> &indexes, std::vector<double> &dists);

		void Search(const Mesh &query_mesh, std::vector<int> &indexes, std::vector<double> &dists);

		/*
		the result vector indexex and dists is depend on the number of the queries and k
		the size of indexes and dists is n * k
		*/
		void Search_knn(const std:: vector<Kd_tree_point> &queries, std::vector<int> &indexes, std::vector<double> &dists, int k);

		void Search_knn(const Mesh &query_mesh, std::vector<int> &indexed, std::vector<double> &dists, int k);

	private:
		CUDA_KDNode *m_gpu_nodes;
		int *m_gpu_indexes;
		Kd_tree_point *m_gpu_points;

		int m_num_points;
	};

	void CheckCUDAError(const char *msg);
}
#endif /* H_CUDA_KDTREE_H */