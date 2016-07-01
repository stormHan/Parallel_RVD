/*

	Parallel compute the Restricted Voronoi Diagram

*/

#include "Command_line.h"
#include "Mesh_io.h"
#include "Mesh.h"
#include "Points.h"
#include "Mesh_repair.h"

#include "RVD.h"

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

	if (!mesh_load(mesh_filename, M_in))
	{
		fprintf(stderr, "cannot load Mesh into M_in!");
		return -1;
	}

	if (!mesh_load(points_filename, points_in, false))
	{
		fprintf(stderr, "cannot load points into points_in!");
		return -1;
	}

	mesh_repair(M_in);
	Points points(points_in);

	/*
		ÉèÖÃKdtree
	*/
	RestrictedVoronoiDiagram *m_RVD = new RestrictedVoronoiDiagram(&M_in, &points);
	m_RVD->compute_RVD();

	return 0;
}
