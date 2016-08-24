/*
	compute the Restricted Voronoi Diagram 
	and Update the seeds

	with Cuda to accelerate the program

	starting time : 2016.7.25
*/

#include "Command_line.h"


using namespace Parallel_RVD;

int main(int argc, char** argv)
{
	std::vector<std::string> filenames;
	if (!Cmd::parse_argu(argc, argv, filenames))
	{
		fprintf(stderr, "cannot parse the argument into filename!");
		return -1;
	}

	std::string mesh_filename = filenames[0];
	std::string points_filename = filenames[0];
	std::string output_filename;
	
	if (filenames.size() >= 2){
		points_filename = filenames[1];
	}
	output_filename = (filenames.size() == 3) ? filenames[2] : "out.eobj";

	

	return 0;
}
