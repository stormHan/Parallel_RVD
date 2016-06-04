#include "Mesh_io.h"

namespace P_RVD
{
	bool mesh_load(const std::string _filepath, Mesh& _M, FileType _filetype)
	{
		
		FILE* file = fopen(_filepath.c_str(), "r");
		if (file == NULL)
		{
			fprintf(stderr, "Impossible to open the file!");
			return false;
		}

		while (true)
		{
			char lineHeader[128];	//we suppose that the length of each line is below 128
			int res = fscanf(file, "%s", lineHeader);
			if (res == EOF)
				break;	//end of the file

			if (strcmp(lineHeader, "v") == 0)
			{
				Vector3f vertex;
				fscanf(file, "%f %f %f", &vertex.x, &vertex.y, &vertex.z);
				_M.meshVertices.addPoint(vertex);
			}
		}

		return true;
	}
}