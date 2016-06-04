/*

	load / save the mesh attributes from file.

*/

#ifndef H_MESH_IO_H
#define H_MESH_IO_H

#include "Common.h"
#include "Mesh.h"


#include <string>

namespace P_RVD
{
	enum FileType
	{
		OBJfile = 0,
		PTXfile = 1
	};

	bool mesh_load(const std::string _filepath, Mesh& _M, FileType _filetype = OBJfile);
}

#endif /* H_MESH_IO_H */