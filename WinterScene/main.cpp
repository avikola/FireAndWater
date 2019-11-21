// Headers
#include <GL/glut.h>
#include <stdio.h>
#include "FluidSolver.h"

// Fluid Solver
FluidSolver *solver;

// Window Width/Height
int width = 800;
int height = 800;


//
// idle
//
void idle()
{
	glutPostRedisplay();
}

//
// display
//
void display()
{
	glutSwapBuffers();
}

//
// main
//
int main(int argc, char** argv)
{
	// Initialize GLUT
	glutInit(&argc, argv);

	// Request double buffer, depth, and color
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);

	// Set Window Size and Position
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(width, height);

	// Create Window
	glutCreateWindow("Winter Wonderland");

	// GLUT Callbacks
	glutIdleFunc(idle);
	glutDisplayFunc(display);

	glutMainLoop();

	return 0;
}