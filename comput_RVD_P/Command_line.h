/*
handle the command argument

*/

#include <vector>
#include <string>

#ifndef H_COMMAND_LINE_H
#define H_COMMAND_LINE_H

namespace P_RVD
{
	namespace Cmd
	{
		/*
		parse the argument.
		take the infomation about the path of input into the filenames

		return
		true : successfully parsed
		*/
		bool parse_argu(int argument_nb, char** argument, std::vector<std::string>& filenames);
	}
}

#endif /* H_COMMAND_LINE_H */