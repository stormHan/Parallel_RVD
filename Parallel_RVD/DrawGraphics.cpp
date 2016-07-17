
#include "DrawGraphics.h"

namespace P_RVD
{
	/*
		set the initial attributes
	*/
	void GraphicsDrawer::Init()
	{
		cameraPosition = { 0.0f, 0.0f, 0.0f };
		centerPosition = { 0.0f, 0.0f, 3.0f };
		upDirection = { 0.0f, 1.0f, 0.0f };

		// set the position of light
		lightPosition[0] = 0.0f;
		lightPosition[1] = 0.0f;
		lightPosition[2] = 3.0f;
		lightPosition[3] = 0.0f;

		//set the ambient, diffuse, specular light attributes
		ambient[0] = 0.2f;
		ambient[1] = 0.2f;
		ambient[2] = 0.2f;
		ambient[3] = 1.0f;

		diffuse[0] = 1.0f;
		diffuse[1] = 1.0f;
		diffuse[2] = 0.0f;
		diffuse[3] = 1.0f;

		specular[0] = 0.0f;
		specular[1] = 0.0f;
		specular[2] = 0.0f;
		specular[3] = 1.0f;
		
		shininess = 40.0f;

		glClearColor(0.1f, 0.2f, 0.1f, 1.0f);
		glEnable(GL_DEPTH_TEST);
		//glDepthFunc(GL_LEQUAL);

		GLfloat white_light[] = { 1.0, 1.0, 1.0, 1.0 };
		GLfloat spot_direction[] = { 0.0, 0.0, -3.0};
		glLightf(GL_LIGHT0, GL_SPOT_CUTOFF, 60.f);
		GLfloat light_model_amb[] = { 0.2, 0.2, 0.2, 1.0 };
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, light_model_amb);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);
		glLightfv(GL_LIGHT0, GL_SPECULAR, white_light);
		glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, spot_direction);
		glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
		glLightf(GL_LIGHT0, GL_SPOT_EXPONENT, 2.0);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);

		glEnable(GL_COLOR_MATERIAL);
//		glShadeModel(GL_SMOOTH);

		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

		glMaterialfv(GL_FLOAT, GL_AMBIENT, ambient);
		glMaterialfv(GL_FLOAT, GL_DIFFUSE, diffuse);
		glMaterialfv(GL_FLOAT, GL_SPECULAR, specular);
		glMaterialf(GL_FLOAT, GL_SHININESS, shininess);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluLookAt(cameraPosition.x, cameraPosition.y, cameraPosition.z,
			centerPosition.x, centerPosition.y, centerPosition.z,
			upDirection.x, upDirection.y, upDirection.z);

		glMatrixMode(GL_MODELVIEW);
	}

	void GraphicsDrawer::Run()
	{
		glFrontFace(GL_CW);
		glCullFace(GL_BACK);
		glEnable(GL_CULL_FACE);

		
	}

	void GraphicsDrawer::DrawPoints(Points _points)
	{
		t_index number = (t_index)_points.getPointsNumber();

		for (t_index i = 0; i < number; ++i)
		{
			glPushMatrix();
			glPointSize(1.0f);
			glBegin(GL_POINTS);
				glVertex3d(_points.getPoint(i).x, _points.getPoint(i).y, _points.getPoint(i).z);
			glEnd();
			glPopMatrix();
		}
	}

	void GraphicsDrawer::DrawMesh(const Mesh& _m)
	{
		t_index facet_nb = _m.meshFacets.getFacetsNumber();
		
		//test
		glPushMatrix();
		glPointSize(10.0f);
		glBegin(GL_POINTS);
		glVertex3d(0.0f, 0.0f, 0.0f);
		glEnd();
		glPopMatrix();
		
		for (t_index i = 0; i < facet_nb; ++i)
		{
			Facet temp_facet = _m.meshFacets.getFacet(i);
			Vector3d pos1 = _m.meshVertices.getPoint(temp_facet.m_v1);
			Vector3d pos2 = _m.meshVertices.getPoint(temp_facet.m_v2);
			Vector3d pos3 = _m.meshVertices.getPoint(temp_facet.m_v3);

			DrawTriangle(pos1, pos2, pos3);
		}
	}

	void GraphicsDrawer::DrawTriangle(Vector3d _v1, Vector3d _v2, Vector3d _v3)
	{
		Vector3d normal = Math::computeNormal(_v1, _v2, _v3);

		glPushMatrix();
		glBegin(GL_TRIANGLES);
			glNormal3d(normal.x, normal.y, normal.z);
			glVertex3d(_v1.x, _v1.y, _v1.z);
			glNormal3d(normal.x, normal.y, normal.z);
			glVertex3d(_v2.x, _v2.y, _v2.z);
			glNormal3d(normal.x, normal.y, normal.z);
			glVertex3d(_v3.x, _v3.y, _v3.z);
		glEnd();
		glPopMatrix();
	}

}