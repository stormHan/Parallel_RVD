
/*
	A class responsible for the Drawing about Geometry
*/

#ifndef H_DRAWGRAPHICS_H
#define H_DRAWGRAPHICS_H

#include <gl\freeglut.h>

#include "Mesh.h"
#include "Points.h"
#include "math_3d.h"

#include "Math_basics.h"

namespace P_RVD
{
	class GraphicsDrawer
	{
	public:
		/*
			Initialize the environmnet of the Grawing
		*/
		void Init();

		/*
			Run the Renderring
		*/
		void Run();

		/*
			Graw Points.
			Include :
			1, Seeds
			2, Intersection Vertex
		*/
		void DrawVertex(Vector3d _point) {};

		/*
			Graw Points
		*/
		void DrawPoints(Points _points);

		/*
			Graw Triangle
			Include :
			1, Mesh facets
			2, RVD clipped Polygon (for each triangle)
		*/
		void DrawMesh(const Mesh& _m);
		
		/*
			Draw Triangle with 3 exact position of vertex
		*/
		void DrawTriangle(Vector3d _v1, Vector3d _v2, Vector3d _v3);
		
	protected:
		
		/*
			camera attributes
		*/
		Vector3f cameraPosition, centerPosition, upDirection;

		/*
			light attributes
		*/
		GLfloat lightPosition[4];

		GLfloat shininess;
		GLfloat ambient[4], diffuse[4], specular[4];

	};
	
}

#endif /* H_DRAWGRAPHICS_H */