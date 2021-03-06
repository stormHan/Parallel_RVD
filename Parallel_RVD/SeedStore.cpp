
#include "SeedStore.h"

namespace P_RVD
{
	void SeedStore::addInformation(SeedWeightPosition _swp, t_index _t)
	{
		seed_updating_map.insert(std::make_pair(_swp, _t));
	
		update_information_nb++;

		if (seed_updating_map.size() != update_information_nb)
		{
			printf("  ");
		}
	}

	void SeedStore::addInformation(Vector3d _pos, double _weight, t_index _t)
	{
		addInformation(SeedWeightPosition(_pos, _weight), _t);
	}

	void SeedStore::UpdateSeeds()
	{

		for (std::map<SeedWeightPosition, t_index>::const_iterator
			iter = seed_updating_map.begin();
			iter != seed_updating_map.end();
		++iter)
		{
			SeedWeightPosition temp = m_seed_weight_pos[iter->second];
			temp.weight += iter->first.weight;
			temp.center += (iter->first.center * iter->first.weight);
			m_seed_weight_pos[iter->second] = temp;
		}
		std::ofstream a("point0_update_info");
		for (std::map<SeedWeightPosition, t_index>::const_iterator
			iter = seed_updating_map.begin();
			iter != seed_updating_map.end();
		++iter)
		{
			if (iter->second == 0)
			{
				a << "x : " << iter->first.center.x << " y: " << iter->first.center.y << " z : " << iter->first.center.z
					<< " w : " << iter->first.weight << std::endl;
			}
		}
		for (std::vector<SeedWeightPosition>::iterator
			iter = m_seed_weight_pos.begin();
			iter != m_seed_weight_pos.end();
		++iter)
		{
			iter->center = iter->center / iter->weight;
		}
	}

}