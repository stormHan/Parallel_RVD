#include "Mesh_io.h"

namespace P_RVD
{
	bool mesh_load(const std::string _filepath, Mesh& _M, bool _meshpoints, FileType _filetype)
	{
		bool hasVt = false, hasVn = false;

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
				Vector3d vertex;
				fscanf(file, "%lf %lf %lf", &vertex.x, &vertex.y, &vertex.z);
				_M.meshVertices.addPoint(vertex);
			}

			if (strcmp(lineHeader, "vt") == 0 && hasVt != true)
				hasVt = true;

			if (strcmp(lineHeader, "vn") == 0 && hasVn != true)
				hasVn = true;

			if (strcmp(lineHeader, "f") == 0 && _meshpoints)
			{
				Vector3i t_vIndex, t_vtIndex, t_vnIndex;
				if (hasVn && hasVt)
				{
					int matched = fscanf(file, "%d/%d/%d %d/%d/%d %d/%d/%d", &t_vIndex.x, &t_vtIndex.x, &t_vnIndex.x,
						&t_vIndex.y, &t_vtIndex.y, &t_vnIndex.y, &t_vIndex.z, &t_vtIndex.z, &t_vnIndex.z);
					if (matched != 9)
					{
						printf("File can not be read by the parser");
						return false;
					}
				}
				else if (hasVn || hasVt)
				{
					int matched = fscanf(file, "%d/%d %d/%d %d/%d", &t_vIndex.x, &t_vtIndex.x, &t_vIndex.y, &t_vtIndex.y,
						&t_vIndex.z, &t_vtIndex.z);
					if (matched != 6)
					{
						printf("File can not be read by the parser");
						return false;
					}
				}
				else
				{
					int matched = fscanf(file, "%d %d %d", &t_vIndex.x, &t_vIndex.y, &t_vIndex.z);
				}
				_M.meshFacets.addFacet((t_index)t_vIndex.x - 1, (t_index)t_vIndex.y - 1, (t_index)t_vIndex.z - 1);
			}
		}

		return true;
	}

	
}