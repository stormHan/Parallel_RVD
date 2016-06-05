/*

	Parallel compute the Restricted Voronoi Diagram


*/

#include "Command_line.h"
#include "Mesh_io.h"
#include "Mesh.h"


int main(int argc, char** argv)
{
	using namespace P_RVD;

	std::vector<std::string> filenames;
	if (!Cmd::parse_argu(argc, argv, filenames))
	{
		fprintf(stderr, "cannot parse the argument into filenames!");
		return -1;
	}
	
	std::string mesh_filename = filenames[0];
	std::string points_filename = filenames[0];
	std::string output_filename;
	if (filenames.size() >= 2) {
		points_filename = filenames[1];
	}
	output_filename = (filenames.size() == 3) ? filenames[2] : "out.eobj";

	Mesh M_in, M_out, points_in;
	FileType mesh_type = OBJfile;
	
	if (!mesh_load(filenames[0], M_in))
	{
		fprintf(stderr, "cannot load Mesh!");
		return -1;
	}

	return 0;
}
