#include "Kdtree.h"
#include <cstdio>
#include <algorithm>
#include <numeric>
#include <cmath>

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