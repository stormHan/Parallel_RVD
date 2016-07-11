/*
	this class is used to store the information to each seed's updating 
	data based on the RVD it(its cell) clipped. 
*/

#ifndef H_SEEDSTORE_H
#define H_SEEDSTORE_H

#include <map>
#include <string>
#include <vector>

#include "math_3d.h"
#include "Common.h"


namespace P_RVD{

	/*
	a data structure to record  the Weight and Position
	of a point.
	*/
	class SeedWeightPosition
	{
	public:
		SeedWeightPosition() : weight(0.00), center(Vector3d(0.0, 0.0, 0.0)) {};

		SeedWeightPosition(Vector3d _pos, double _weight) : weight(_weight), center(_pos) {};

		double weight;
		Vector3d center;

		bool operator < (const SeedWeightPosition _swp) const
		{
			if (center.x < _swp.center.x) 
				return true;
			else if (center.x == _swp.center.x && center.y < _swp.center.y)
			{
				return true;
			}
			else if (center.x == _swp.center.x && center.y == _swp.center.y && center.z < _swp.center.z)
			{
				return true;
			}
			else
				return false;
		}
	};

	class SeedStore
	{
	public:
		SeedStore() : update_information_nb(0) {};

		SeedStore(t_index _seedsnb) : update_information_nb(0), seed_nb(_seedsnb) {};

		/*
			set the seed_weight_pos vector
			must be done before using the object
		*/
		void setPositionVector()
		{
			if (seed_nb == 0)
			{
				fprintf(stderr, "error setting the position vector!/n");
				exit(0);
			}
			m_seed_weight_pos = std::vector<SeedWeightPosition>(seed_nb, SeedWeightPosition(
				Vector3d(0.0, 0.0, 0.0), 0));
		}

		/*
			set seeds number
		*/ 
		void setSeedsNumber(t_index _t) { seed_nb = _t; };

		/*
		clear the information of SeedStore
		*/
		void clear()
		{
			seed_updating_map.clear();
			update_information_nb = 0;
		}

		/*
			push new pair to the map

			_swp ; the weight and position information to a seed
			_t   : the indice of a seed.
		*/
		void addInformation(SeedWeightPosition _swp, t_index _t);

		void addInformation(Vector3d _pos, double _weight, t_index _t);

		/*
			get the number of the updateing information
		*/
		t_index getSeedUpdatingInformation(){ return update_information_nb; };

		void UpdateSeeds();

	protected:
		/*
			a map to record the Polygon information (SeedWightPosition) hash
			the exact seed.
			*/
		std::map<SeedWeightPosition, t_index> seed_updating_map;
		t_index update_information_nb;

		t_index seed_nb;
		std::vector<SeedWeightPosition> m_seed_weight_pos;
	};
}

#endif /* H_SEEDSTORE_H */