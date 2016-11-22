#include "AnnKdtree.h"

namespace P_RVD{
	m_AnnKdtree::m_AnnKdtree(const double* data, int n){
		m_eps = 0;
		m_points = annAllocPts(n, 3);
		memcpy(&m_points[0][0], data, sizeof(double) * 3 * n);
		m_kdtree = new ANNkd_tree(m_points, n, 3);
	}

	m_AnnKdtree::~m_AnnKdtree(){
		annDeallocPts(m_points);
		delete m_kdtree;
		annClose();
	}

	void m_AnnKdtree::queryNearestNeighbors(const double pt[3], int k, std::vector<int>& indexes)
	{
		std::vector<ANNdist> m_dist(k);
		indexes.resize(k);
		double xyz[3];
		xyz[0] = pt[0]; xyz[1] = pt[1]; xyz[2] = pt[2];
		m_kdtree->annkSearch(xyz, k, &indexes[0], &m_dist[0], m_eps);
	}

	void m_AnnKdtree::queryNearestNeighbors(Vector3d pt, int k, std::vector<int>& indexes)
	{
		std::vector<ANNdist> m_dist(k);
		indexes.resize(k);
		double xyz[3];
		xyz[0] = pt.x; xyz[1] = pt.y; xyz[2] = pt.z;
		m_kdtree->annkSearch(xyz, k, &indexes[0], &m_dist[0], m_eps);
	}

}