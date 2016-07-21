
/*
	Create the Windows
	with GLUT lib
*/

#ifndef H_GLUT_GENERATOR_H
#define H_GLUT_GENERATOR_H

#include <GL\freeglut.h>

#include "Common.h"

namespace P_RVD
{
	static bool sWithDepth = false;
	static bool sWithStencil = false;

	void GLUTBackendInit(int argc, char** argv, bool WithDepth, bool WithStencil)
	{
		sWithDepth = WithDepth;
		sWithStencil = WithStencil;

		glutInit(&argc, argv);

		uint DisplayMode = GLUT_DOUBLE | GLUT_RGBA;

		if (WithDepth) {
			DisplayMode |= GLUT_DEPTH;
		}

		if (WithStencil) {
			DisplayMode |= GLUT_STENCIL;
		}

		glutInitDisplayMode(DisplayMode);

		glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
	}

	bool GLUTBackendCreateWindow(unsigned int Width, unsigned int Height, bool isFullScreen, const char* pTitle)
	{
		if (isFullScreen) {
			char ModeString[64] = { 0 };
			int bpp = 32;
			_snprintf_s(ModeString, sizeof(ModeString), "%dx%d:%d@60", Width, Height, bpp);
			glutGameModeString(ModeString);
			glutEnterGameMode();
		}
		else {
			glutInitWindowPosition(0, 0);
			glutInitWindowSize(Width, Height);
			glutCreateWindow(pTitle);
		}
		return true;
	}

	//static void InitCallbacks()
	//{
	//	glutDisplayFunc(RenderSceneCB);
	//	glutIdleFunc(IdleCB);
	//	glutSpecialFunc(SpecialKeyboardCB);
	//	glutPassiveMotionFunc(PassiveMouseCB);
	//	glutKeyboardFunc(KeyboardCB);
	//	glutMouseFunc(MouseCB);
	//}
}


#endif /* H_GLUT_GENERATOR_H */