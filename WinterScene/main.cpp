// Headers
#include <GL/glut.h>
#include <stdio.h>
#include "FluidSolver.h"

// Fluid Solver
FluidSolver *solver;

// Mouse Position
int initMousePos[2] = { 0, 0 };
int mousePos[2] = { 0, 0 };

// Mouse Button Press States
int leftMouseButton = 0;
int middleMouseButton = 0;
int rightMouseButton = 0;

// Window Width/Height
int width = 800;
int height = 800;

// Velocity Line Length
float velLength = 5.0f;

//
// drawVelocities
//
void drawVelocities()
{
	// Get Position/Velocity Fields
	point* points = solver->getPositions();
	velocity* vels = solver->getVelocities();

	// Draw Style
	glColor3f(0.0f, 1.0f, 0.0f);
	glLineWidth(1.0f);

	// Draw Velocity Lines
	glBegin(GL_LINES);
	for (int i = 0; i < solver->getSize(); ++i)
	{
		// Get Position/Velocity
		point p = points[i];
		velocity v = vels[i];

		// Draw
		glVertex2f(p.x, p.y);
		glVertex2f(p.x + (v.x * velLength), p.y + v.y * velLength);
	}
	glEnd();
}

//
// mouseButton
//
void mouseButton(int button, int state, int x, int y)
{
	// Set Mouse Position
	mousePos[0] = x;
	mousePos[1] = y;

	// Set Initial Mouse Position
	initMousePos[0] = x;
	initMousePos[1] = y;

	// Get Mouse Button State
	switch (button)
	{
		case GLUT_LEFT_BUTTON:
			leftMouseButton = (state == GLUT_DOWN);
			break;
		case GLUT_MIDDLE_BUTTON:
			middleMouseButton = (state == GLUT_DOWN);
			break;
		case GLUT_RIGHT_BUTTON:
			rightMouseButton = (state == GLUT_DOWN);
			break;
	}
}

//
// getMouseInput
//
void getMouseInput()
{
	// Reset the Initial Velocities
	solver->resetInitialVelocities();

	// Get Rows/Cols Size
	int rows = solver->getRows();
	int cols = solver->getCols();

	// Position in Grid
	int xPos; int yPos;

	// If Right Mouse Button Pressed
	if (rightMouseButton)
	{
		// Get Position relative to Grid
		xPos = (int)((float)(mousePos[0]) / width * (rows));
		yPos = (int)((float)(height - mousePos[1]) / height * (cols));

		// Check if within bounds of Grid
		if ((xPos > 0) && (xPos < rows - 1) && (yPos > 0) && (yPos < cols - 1))
		{
			// If Right Mouse Button Pressed
			if (rightMouseButton)
			{
				// Calculate velocities
				float xVel = 1.0f * (mousePos[0] - initMousePos[0]);
				float yVel = 1.0f * (initMousePos[1] - mousePos[1]);

				// Set initial velocity
				solver->setInitVelocity(xPos, yPos, xVel, yVel);
			}

			// Update Initial Mouse Position
			initMousePos[0] = mousePos[0];
			initMousePos[1] = mousePos[1];
		}

		solver->addSource();
	}
}

/**
 * mouseDrag
 */
void mouseDrag(int x, int y)
{
	// Update Mouse Position
	mousePos[0] = x;
	mousePos[1] = y;
}

//
// idle
//
void idle()
{
	// Make the screen update
	glutPostRedisplay();
}

//
// display
//
void display()
{
	// Get Mouse Input
	getMouseInput();

	// Step the Velocity
	solver->stepVelocity(0.0, 1.0);

	// Set Background Color
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// Draw the Velocity
	drawVelocities();

	// Style of Grid
	glColor3f(1.0f, 0.0f, 0.0f);
	glPointSize(1.0f);

	// Get Points
	point* points = solver->getPositions();

	// Double buffer flush
	glutSwapBuffers();
}

/*
 * reshape - Called every time window is resized
 *           to update the projection matrix, and
 *           to preserve aspect ratio
 */
void reshape(int w, int h)
{
	// Initialize Aspect
	GLfloat aspect = (GLfloat)w / (GLfloat)h;

	// Setup image size
	glViewport(0, 0, w, h);

	// Set Perspective
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0f, (float)(solver->getRows()), 0.0f, (float)(solver->getCols()));

	// Set back to ModelView
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

//
// main
//
int main(int argc, char** argv)
{
	// Initialize GLUT
	glutInit(&argc, argv);

	// Initialize Fluid Solver
	solver = new FluidSolver();
	solver->init(800,800);
	solver->reset();

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
	glutMotionFunc(mouseDrag);
	glutMouseFunc(mouseButton);
	glutReshapeFunc(reshape);

	glutMainLoop();

	return 0;
}