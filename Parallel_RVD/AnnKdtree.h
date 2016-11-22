#include <ANN\ANN.h>
#include "math_3d.h"
#include <vector>

namespace P_RVD{
	class m_AnnKdtree{

	public:
		m_AnnKdtree(const double* data, int n);
		~m_AnnKdtree();

		void queryNearestNeighbors(const double pt[3], int k, std::vector<int>& indexes);
		void queryNearestNeighbors(Vector3d pt, int k, std::vector<int>& indexes);

	private:
		ANNkd_tree*		m_kdtree;
		ANNpointArray	m_points;
		double			m_eps;
	};
}