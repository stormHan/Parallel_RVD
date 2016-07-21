/*

	Parallel compute the Restricted Voronoi Diagram

*/

//RVD computing part
#include "Command_line.h"
#include "Mesh_io.h"
#include "Mesh.h"
#include "Points.h"
#include "Mesh_repair.h"

// Renderring part
#include "Glut_generator.h"
#include "DrawGraphics.h"

#include "RVD.h"

#define WINDOWS_WIDTH 1000
#define WINDOWS_HEIGHT 800

using namespace P_RVD;

GraphicsDrawer* m_GraphicsDrawer = new GraphicsDrawer();
Points p_in, p_out;
Mesh M_in;
void RenderCB(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0.1f, 0.1f, 0.2f, 0.0f);

	
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_LIGHTING);
	glColor3f(1.0f, 0.0f, 0.0f);
	m_GraphicsDrawer->DrawPoints(p_in);

	glColor3f(0.0f, 1.0f, 0.0f);
	m_GraphicsDrawer->DrawPoints(p_out);

	glEnable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);
	glColor3f(0.8f, 0.8f, 0.8f);
	glPolygonMode(GL_FRONT, GL_LINE);
	m_GraphicsDrawer->DrawMesh(M_in);

	glFlush();
	glutSwapBuffers();
	glutPostRedisplay();
}

void myResize(int width, int height)
{

	//get the viewport of GLUI and then set the viewport back to that after the resize
	glViewport(0, 0, width, height);//viewport函数用于在窗口中设定一个视角的大小，可以用来将一个窗口划分为多个视角
	glMatrixMode(GL_PROJECTION);

	glPushMatrix();
	glLoadIdentity();
	gluPerspective(60.0, (float)width / height, 0.5f, 20.0f);//重新设定视角
	glPopMatrix();
	gluLookAt(1.0f, 0.0f, 0.0f, 0.0f, 0.f, 0.0f, 0.0f, 1.0f, 0.0f);
	glClearColor((GLclampf)1 / 255, (GLclampf)1 / 255, (GLclampf)1 / 255, 0.0);//将背景刷成灰色
	glEnable(GL_DEPTH_TEST);
}

int main(int argc, char** argv)
{
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

	Mesh  M_out, points_in;
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
	Points points_out(points_in);

	/*
	设置Kdtree
	*/
	RestrictedVoronoiDiagram *m_RVD = new RestrictedVoronoiDiagram(&M_in, &points_out);
	m_RVD->compute_RVD();

	/*
		Set the Windows
		Use the GLUT lib
	*/
	GLUTBackendInit(argc, argv, true, false);
	GLUTBackendCreateWindow(WINDOWS_WIDTH, WINDOWS_HEIGHT, false, "Rvd");

	//GraphicsDrawer* m_GraphicsDrawer = new GraphicsDrawer();
	m_GraphicsDrawer->Init();
	m_GraphicsDrawer->Run();

	p_in = points;
	p_out = points_out;

	glutDisplayFunc(&RenderCB);
	//glutReshapeFunc(&myResize);

	glutPostRedisplay();
	glutMainLoop();

	delete m_GraphicsDrawer;

	return 0;
}
