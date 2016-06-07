/*
	Mesh_repair file

	as the input mesh may has some unexpected characteristics
	we should optimate the mesh
*/

#ifndef H_MESH_REPAIR_H
#define H_MESH_REPAIR_H

#include "Common.h"

namespace P_RVD
{
	class Mesh;

	bool mesh_repair(Mesh& _M);
}

#endif /* H_MESH_REPAIR_H */