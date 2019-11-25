// System Libraries
#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <windows.h>

// OpenGL Libraries
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

// Headers
#include "particle_system.h"
//#include "SOIL.h"
#include "constants.h"
#include "cloud.h"


// Namespaces
using namespace std;

particle_system p(NUMBER_OF_PARTICLES);
float winWidth = 400, winHeight = 400;

//Called when a key is pressed
void handle_keypress(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'a':
	case 'A':
		p.add_particles(10);
		break;

	case 'd':
	case 'D':
		p.delete_particles(10);
		break;


	case 27:
		exit(0);
		break;
	}
	glutPostRedisplay();
}



//Initializes 3D rendering
void init()
{
	// Make big points and wide lines
	glPointSize(1);

	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

	//Enable transparency
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}


//Called when the window is resized
void handle_resize(int w, int h)
{
	//setup windows width and height
	winWidth = (w == 0) ? 1 : w;
	winHeight = (h == 0) ? 1 : h;

	//Tell OpenGL how to convert from coordinates to pixel values
	glViewport(0, 0, winWidth, winHeight);

	//Switch to setting the camera perspective
	glMatrixMode(GL_PROJECTION);

	//Set the camera perspective
	glLoadIdentity(); //Reset the camera
	glOrtho(-LENGTH, LENGTH, -LENGTH, LENGTH, -LENGTH, LENGTH);

	glMatrixMode(GL_MODELVIEW);
}


//Draws the 3D scene
void draw()
{
	glLoadIdentity();

	//Draw particles
	glPushMatrix();
	p.advance(DELTA);
	p.delete_particle();
	p.draw();
	glPopMatrix();

	//Draw overlaying quad for trail
	glColor4f(0, 0, 0, 0.1);
	glBegin(GL_QUADS);
	glVertex3f(-LENGTH, -LENGTH, 100);
	glVertex3f(LENGTH, -LENGTH, 100);
	glVertex3f(LENGTH, LENGTH, 100);
	glVertex3f(-LENGTH, LENGTH, 100);
	glEnd();

	glutSwapBuffers();
	glutPostRedisplay();
}


//Handle mouse movement
void mouse_movement(int x, int y) {
	float ww_ratio = float(x) / winWidth;
	float wh_ratio = float(y) / winHeight;

	p.set_gravity(vec3d((2 * ww_ratio - 1)*LENGTH, (1 - 2 * wh_ratio)*LENGTH, 0));
}


//
// main
//
int main(int argc, char** argv)
{
	srand(time(0));
	p.set_gravity();

	//Initialize GLUT
	glutInit(&argc, argv);

	// Request double buffer, and color
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

	// Set Window Size and Position
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(winWidth, winHeight);

	//Create the window and initialize OpenGL
	glutCreateWindow("Realistic Fire");
	init();

	//Set handler functions for drawing, keypresses, and window resizes
	glutDisplayFunc(draw);
	glutKeyboardFunc(handle_keypress);
	glutReshapeFunc(handle_resize);
	glutPassiveMotionFunc(mouse_movement);


	glutMainLoop();
	return 0; //This line is never reached
}