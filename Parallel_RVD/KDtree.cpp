#include "Kdtree.h"
#include <cstdio>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace P_RVD{
	// Eww global. Need this for the sort function
	// A pointer to the class itself, allows the sort function to determine the splitting axis to sort by
	static KDtree *myself = NULL;

	KDtree::KDtree()
	{
		myself = this;
		m_id = 0;
	}

	KDtree::~KDtree()
	{
		// Delete all the ndoes
		std::vector <KDNode*> to_visit;

		to_visit.push_back(m_root);

		while (to_visit.size()) {
			std::vector <KDNode*> next_search;

			while (to_visit.size()) {
				KDNode *cur = to_visit.back();
				to_visit.pop_back();

				if (cur->left)
					next_search.push_back(cur->left);

				if (cur->right)
					next_search.push_back(cur->right);

				delete cur;
			}

			to_visit = next_search;
		}

		m_root = NULL;
	}

	void KDtree::Create_kdtree(std::vector<Kd_tree_point> &pts, int max_levels)
	{
		m_pts = &pts;
		m_levels = max_levels;

		m_root = new KDNode();
		m_root->id = m_id++;
		m_root->index.resize(pts.size());

		//init the index
		for (unsigned int i = 0; i < pts.size(); ++i)
		{
			m_root->index[i] = i;
		}

		std::vector<KDNode*> to_visit;
		to_visit.push_back(m_root);

		while (to_visit.size())
		{
			std::vector<KDNode*> next_search;
			while (to_visit.size())
			{
				KDNode *node = to_visit.back();
				to_visit.pop_back();;

				if (node->level < max_levels)
				{
					//if current node's index size > 1
					//it has more than 1 point in this space
					//So split it
					if (node->index.size() > 1)
					{
						KDNode *left = new KDNode();
						KDNode *right = new KDNode();

						left->id = m_id++;
						right->id = m_id++;

						Split(node, left, right);

						//Claer current index
						{
							//std::vector<int> dummy;
							//node->index.swap(dummy);
						}

						node->left = left;
						node->right = right;

						node->_left = left->id;
						node->_right = right->id;

						if (left->index.size())
							next_search.push_back(left);
						
						if (right->index.size())
							next_search.push_back(right);
					}
				}
			}
			to_visit = next_search;
		}
	}

	/*
		compare the data of 2 points with the index
	*/
	bool KDtree::SortPoints(const int a, const int b)
	{
		std::vector<Kd_tree_point> &pts = *myself->m_pts;

		return pts[a].coords[myself->m_current_axis] < pts[b].coords[myself->m_current_axis];
	}

	/*
		Split the current node into 2 part: left and right

		we assume that the left and right nodes are already created
	*/
	void KDtree::Split(KDNode *cur, KDNode *left, KDNode *right)
	{
		std::vector<Kd_tree_point> &pts = *m_pts;
		m_current_axis = cur->level % KDTREE_DIM;

		//sort the points by the axis value
		std::sort(cur->index.begin(), cur->index.end(), KDtree::SortPoints);

		int mid = cur->index[cur->index.size() / 2];
		cur->split_value = pts[mid].coords[m_current_axis];

		left->parent = cur;
		right->parent = cur;

		left->level = cur->level + 1;
		right->level = cur->level + 1;

		left->_parent = cur->id;
		right->_parent = cur->id;

		for (unsigned int i = 0; i < cur->index.size(); ++i)
		{
			int idx = cur->index[i];
			if (pts[idx].coords[m_current_axis] < cur->split_value)
				left->index.push_back(idx);
			else
				right->index.push_back(idx);
		}
	}

	/*
		search the best node
		return the index and distance square of the closed point in this node
	*/
	void KDtree::SearchAtNode(KDNode *cur, const Kd_tree_point &query, int *ret_index, double *ret_dist, KDNode **ret_node)
	{
		int best_idx = 0;
		// if the range the coords is larger than 100,best_dist should be adjusted
		double best_dist_Squa = 100.0;

		std::vector<Kd_tree_point> &pts = *m_pts;

		//First pass
		while (true)
		{
			int split_axis = cur->level % KDTREE_DIM;
			m_cmps++;

			if (cur->left == NULL)
			{
				*ret_node = cur;
				
				for (unsigned int i = 0; i < cur->index.size(); ++i)
				{
					m_cmps++;

					int idx = cur->index[i];

					double distSqua = DistanceSquare(query, pts[idx]);
					if (distSqua < best_dist_Squa)
					{
						best_dist_Squa = distSqua;
						best_idx = idx;
					}
				}
				break;
			}
			else if (query.coords[split_axis] < cur->split_value)
			{
				cur = cur->left;
			}
			else
			{
				cur = cur->right;
			}
		}
		*ret_index = best_idx;
		*ret_dist = best_dist_Squa;
	}

	void KDtree::SearchAtNodeRange(KDNode *cur, const Kd_tree_point &query, double range, int *ret_index, double *ret_dist)
	{
		int best_idx = 0;
		double best_dist = 100.0;

		std::vector<Kd_tree_point> &pts = *m_pts;
		std::vector<KDNode* > to_visit;

		to_visit.push_back(cur);
		while (to_visit.size() > 0)
		{
			std::vector<KDNode*> next_search;
			while (to_visit.size())
			{
				cur = to_visit.back();
				to_visit.pop_back();

				int split_axis = cur->level % KDTREE_DIM;

				if (cur->left == NULL)
				{
					for (unsigned int i = 0; i < cur->index.size(); i++) {
						m_cmps++;

						int idx = cur->index[i];
						double d = DistanceSquare(query, pts[idx]);

						if (d < best_dist) {
							best_dist = d;
							best_idx = idx;
						}
					}
				}
				else
				{
					double d = query.coords[split_axis] - cur->split_value;

					// There are 3 possible scenarios
					// The hypercircle only intersects the left region
					// The hypercircle only intersects the right region
					// The hypercricle intersects both

					m_cmps++;

					if (fabs(d) > range) {
						if (d < 0)
							next_search.push_back(cur->left);
						else
							next_search.push_back(cur->right);
					}
					else {
						next_search.push_back(cur->left);
						next_search.push_back(cur->right);
					}
				}
			}
		}

	}

	void KDtree::Search(const Kd_tree_point &query, int *ret_index, double *ret_dist)
	{
		//First pass
		//Find the first closest node, this will be the upper bound for the next searches
		std::vector<Kd_tree_point> &pts = *m_pts;
		KDNode *best_node = NULL;
		int best_idx = 0;
		double best_dist = 100.0;
		double radius = 0.0;
		m_cmps = 0;

		SearchAtNode(m_root, query, &best_idx, &best_dist, &best_node);
		radius = sqrt(best_dist);

		//Second pss
		//Find ohter possible candidates
		KDNode *cur = best_node;
		while (cur->parent != NULL)
		{
			//Go up
			KDNode *parent = cur->parent;
			int split_axis = (parent->level) % KDTREE_DIM;

			//Search the other node
			int temp_idx;
			double temp_dist = 100.0;
			KDNode *temp_node;
			KDNode *search_node = NULL;

			if (fabs(parent->split_value - query.coords[split_axis]) <= radius)
			{
				//Search opposite node
				if (parent->left != cur)
					SearchAtNodeRange(parent->left, query, radius, &temp_idx, &temp_dist);
				else
					SearchAtNodeRange(parent->right, query, radius, &temp_idx, &temp_dist);
			}

			if (temp_dist < best_dist)
			{
				best_idx = temp_idx;
				best_dist = temp_dist;
			}
			cur = parent;
		}
		*ret_index = best_idx;
		*ret_dist = best_dist;
	}

	inline double KDtree::DistanceSquare(const Kd_tree_point &a, const Kd_tree_point &b) const
	{
		double dist = 0.0;

		for (int i = 0; i < KDTREE_DIM; ++i)
		{
			double d = a.coords[i] - b.coords[i];
			dist += d * d;
		}
		return dist;
	}
}