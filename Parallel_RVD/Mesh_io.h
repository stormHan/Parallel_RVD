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

	/*
		load the mesh attributes from the file

		meshpoints : true	-- load a mesh
					 flase	-- load points, just vertices
	*/
	bool mesh_load(const std::string _filepath, Mesh& _M, bool _meshpoints = true, FileType _filetype = OBJfile);

	/*
		save the mesh.
	*/

	
}

#endif /* H_MESH_IO_H */