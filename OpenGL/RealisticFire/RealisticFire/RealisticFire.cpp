// Standard Libraries
#include <time.h>

// OpenGL Libraries
#include <GL/glut.h>

// Headers
#include "particle_system.h"
#include "Fluid.h"

// Create Fluid Solver
FluidSolver *fluidSolver;

// Create Particle System
particle_system p(NUMBER_OF_PARTICLES);

int disp_type = 1;

int win_x = 800;
int win_y = 800;

int mouse_down[3];
int omx;
int omy;
int mx;
int my;

/**
 * drawVelocity
 */
void drawVelocity()
{	
	// Get Points
	point* pt = fluidSolver->getPoints();

	// Get Velocities
	float *vx = fluidSolver->getVX();
	float *vy = fluidSolver->getVY();

	// Line Color/Style
	glColor3f(0.0f, 0.0f, 1.0f);
	glLineWidth(1.0f);

	// Draw Velocity Lines
	glBegin(GL_LINES);
		for (int i = 0; i < fluidSolver->getTotSize(); i++)
		{
			glVertex2f(pt[i].x, pt[i].y);
			glVertex2f(pt[i].x + vx[i] * 10.0f, pt[i].y + vy[i] * 10.0f);
		}
	glEnd();
}

/**
 * drawDensity
 */
void drawDensity()
{
	float x;
	float y;
	float d00;
	float d01;
	float d10;
	float d11;

	int rowSize = fluidSolver->getRows();
	int colSize = fluidSolver->getColSize();

	glBegin(GL_QUADS);
	for (int i = 1; i <= rowSize - 2; i++)
	{
		x = (float)i;
		for (int j = 1; j <= colSize - 2; j++)
		{
			y = (float)j;

			d00 = fluidSolver->getDens(i, j);
			d01 = fluidSolver->getDens(i, j + 1);
			d10 = fluidSolver->getDens(i + 1, j);
			d11 = fluidSolver->getDens(i + 1, j + 1);

			glColor3f(0.0f + d00, 0.0f + d00, 0.0f + d00); glVertex2f(x, y);
			glColor3f(0.0f + d10, 0.0f + d10, 0.0f + d10); glVertex2f(x + 1.0f, y);
			glColor3f(0.0f + d11, 0.0f + d11, 0.0f + d11); glVertex2f(x + 1.0f, y + 1.0f);
			glColor3f(0.0f + d01, 0.0f + d01, 0.0f + d01); glVertex2f(x, y + 1.0f);
		}
	}
	glEnd();
}

/**
 * drawFire
 */
void drawFire()
{
	glLoadIdentity();

	//Draw particles
	glPushMatrix();
	p.advance(DELTA);
	p.delete_particle();
	p.draw();
	glPopMatrix();
}

/**
 * getMouseInput
 */
void getMouseInput()
{
	// Reset initial velocities/densities
	fluidSolver->resetInitialFields();

	// Get Rows/Cols
	int rowSize = fluidSolver->getRows();
	int colSize = fluidSolver->getColSize();

	// Check if Left/Right Mouse Button Clicked
	if (mouse_down[0] || mouse_down[2])
	{
		// Get Mouse position relative to grid
		int xPos = (int)((float)(omx) / win_x * (rowSize));
		int yPos = (int)((float)(win_y - omy) / win_y * (colSize));

		// Bounds Check
		if (xPos > 0 && xPos < rowSize - 1 && yPos > 0 && yPos < colSize - 1)
		{
			// Left Mouse Button
			if (mouse_down[0])
			{
				// Get the velocity values
				float xVel = 1.0f * (mx - omx);
				float yVel = 1.0f * (omy - my);

				// Set the initial velocity
				fluidSolver->setInitVelocity(xPos, yPos, xVel, yVel);
			}

			// Right Mouse Button
			if (mouse_down[2])
			{
				// Density Value
				float density = 10.0f;
				fluidSolver->setInitDensity(xPos, yPos, density);
			}

			omx = mx;
			omy = my;
		}

		// Set Velocities/Densities
		fluidSolver->addSource();
	}
}

void key_func(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'v':
	case 'V':
		disp_type = (disp_type + 1) % 2;
		break;

	case ' ':
		if (fluidSolver->isRunning() == 1)
		{
			fluidSolver->stop();
		}
		else
		{
			fluidSolver->start();
		}
		break;
	case 'c':
	case 'C':
		fluidSolver->resetFields();
		break;
	case 27: // escape
		exit(0);
		break;
	}
}

void mouse_func(int button, int state, int x, int y)
{
	omx = x;
	omy = y;

	mx = x;
	my = y;

	mouse_down[button] = state == GLUT_DOWN;
}

void motion_func(int x, int y)
{
	mx = x;
	my = y;
}

/**
 * reshape
 */
void reshape(int w, int h)
{
	// Set Image Size
	glViewport(0, 0, win_x, win_y);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// Set Camera Perspective
	glOrtho(0.0f, (float)(fluidSolver->getRows()) - BUFFER, 0.0f, (float)(fluidSolver->getColSize()) - BUFFER, -LENGTH, LENGTH);

	// Set back to Model View
	glMatrixMode(GL_MODELVIEW);
}

/**
 * idle - Callback that is called
 *        when no other events to be
 *        handled
 */
void idle()
{
	// Make the screen update
	glutPostRedisplay();
}

/**
 * display
 */
void display()
{
	// Get input from Mouse
	getMouseInput();

	// Move the Velocity/Density forward 1 timestep
	fluidSolver->stepVelocity();
	fluidSolver->stepDensity();

	// Draw Density/velocity
	if (disp_type == 0) drawDensity();
	if (disp_type == 1) drawVelocity();

	// Draw the Fire
	drawFire();

	// Double Buffer Flush
	glutSwapBuffers();
}


/*
 * init
 */
void init()
{
	// Make big points and wide lines
	glPointSize(1);

	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

	//Enable transparency
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Set background color [Black]
	glClearColor(0.0, 0.0, 0.0, 0.0); 
}

/**
 * main
 */
int main(int argc, char** argv)
{
	// Particle Initialization
	srand(time(0));
	p.set_gravity();

	// Initialize Fluid Solver
	fluidSolver = new FluidSolver(100, 100, 1.0, 0.0f);
	fluidSolver->init();
	fluidSolver->resetFields();

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

	// Set Window size / Position
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(win_x, win_y);

	// Create Window
	glutCreateWindow("Realistic Fire");

	// GLUT Callbacks
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(key_func);
	glutMouseFunc(mouse_func);
	glutMotionFunc(motion_func);

	// Initialize States
	init();

	glutMainLoop();

	return 0;
}