/*
	Kdtree 
	for find the nearest points
*/
#ifndef H_KDTREE_H
#define H_KDTREE_H

#include <vector>

#define KDTREE_DIM 3 //DATA dimension


namespace P_RVD{
	/*
		simple struct to store the data of kdtree
		searching point
		*/
	struct Kd_tree_point
	{
		double coords[3];
	};

	class KDNode
	{
	public:
		KDNode()
		{
			parent = NULL;
			left = NULL;
			right = NULL;

			split_value = -1;
			_parent = -1;
			_right = -1;
			_left = -1;
		}

		int id; // for GPU 
		int level;
		KDNode *parent, *left, *right;
		double split_value;
		int _parent, _left, _right; // for GPU
		std::vector<int> index; //indices to points
	};

	class KDtree
	{
	public:
		KDtree();
		~KDtree();

		/*
			Create the Kree in CPU with the pts
			*/
		void Create_kdtree(std::vector<Kd_tree_point> &pts, int max_levels = 99);

		void Search(const Kd_tree_point &query, int *ret_index, double *ret_sq_dist);
		int GetNumNodes() const { return m_id; }
		KDNode* GetRoot() const { return m_root; }

		static bool SortPoints(const int a, const int b);

	private:
		std::vector <Kd_tree_point> *m_pts;
		KDNode *m_root;
		int m_current_axis;
		int m_levels;
		int m_cmps; // count how many comparisons were made in the tree for a query
		int m_id; // current node ID

		void Split(KDNode *cur, KDNode *left, KDNode *right);
		void SearchAtNode(KDNode *cur, const Kd_tree_point &query, int *ret_index, double *ret_dist, KDNode **ret_node);
		void SearchAtNodeRange(KDNode *cur, const Kd_tree_point &query, double range, int *ret_index, double *ret_dist);
		inline double DistanceSquare(const Kd_tree_point &a, const Kd_tree_point &b) const;
	};
}
#endif /* H_KDTREE_H */