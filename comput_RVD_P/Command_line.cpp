
#include "Command_line.h"

namespace Parallel_RVD
{
	namespace Cmd
	{
		bool parse_argu(int argument_nb, char** argument, std::vector<std::string>& filenames)
		{
			if (argument_nb <= 1)
			{
				fprintf(stderr, "no argument pointed to filename!");
				return false;
			}

			if (filenames.size() >= 1)
			{
				fprintf(stderr, "already have filenames!");
				return false;
			}
			for (int i = 1; i < argument_nb && i <= 3; ++i)
			{
				filenames.push_back(argument[i]);
			}

			return true;
		}
	}
}