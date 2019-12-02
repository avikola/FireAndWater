// Standard Libraries
#include <time.h>

// OpenGL Libraries
#include <GL/glut.h>

// Headers
#include "particle_system.h"
#include "Fluid.h"
#include "SOIL.h"

// Create Fluid Solver
FluidSolver *fluidSolver;

// Window Width/Height
const int WIDTH = 800;
const int HEIGHT = 800;
GLuint tex;

// Create Particle System
particle_system p(NUMBER_OF_PARTICLES);

// Display Mode
int displayMode = 0;

// Current Mouse Position
int _mousePos[2] = { 0, 0 };
int _origMousePos[2] = { 0, 0 };

//  Mouse Buttons
int _leftMouseButton = 0;
int _middleMouseButton = 0;
int _rightMouseButton = 0;

// Initial Velocities
int x_velocity = 0;		// x velocity
int y_velocity = 60;	// y velocity

// X Pos, Y Pos
static int xPos;
static int yPos;

/**
 * drawVelocity
 */
void drawVelocity()
{	
	// Get Points
	point* p = fluidSolver->getPoints();

	// Get Velocities
	velocity* v = fluidSolver->getVelocities();

	// Line Color/Style
	glColor3f(0.0f, 0.0f, 1.0f);
	glLineWidth(1.0f);

	// Draw Velocity Lines
	glBegin(GL_LINES);
		for (int i = 0; i < fluidSolver->getSize(); i++)
		{
			glVertex2f(p[i].x, p[i].y);
			glVertex2f(p[i].x + v[i].x, p[i].y + v[i].y);
		}
	glEnd();
}

/**
 * drawDensity
 */
void drawDensity()
{
	// Get Rows/Cols
	int rowSize = fluidSolver->getRows();
	int colSize = fluidSolver->getCols();

	// Draw Fluid Density
	glBegin(GL_QUADS);
		for (int i = 1; i <= rowSize - 2; i++)
		{
			float x = (float)i;
			for (int j = 1; j <= colSize - 2; j++)
			{
				float y = (float)j;

				float d00 = fluidSolver->getDensity(i, j);
				float d01 = fluidSolver->getDensity(i, j + 1);
				float d10 = fluidSolver->getDensity(i + 1, j);
				float d11 = fluidSolver->getDensity(i + 1, j + 1);

				glColor3f(0.0f + d00, 0.0f + d00, 0.0f + d00); glVertex2f(x, y);
				glColor3f(0.0f + d10, 0.0f + d10, 0.0f + d10); glVertex2f(x + 1.0f, y);
				glColor3f(0.0f + d11, 0.0f + d11, 0.0f + d11); glVertex2f(x + 1.0f, y + 1.0f);
				glColor3f(0.0f + d01, 0.0f + d01, 0.0f + d01); glVertex2f(x, y + 1.0f);
			}
		}
	glEnd();
}


/**
 * Initialize fire texture.
 */
void initTex() {
	int width = 280;
	int height = 280;
	glGenTextures(1, &tex);
	glBindTexture(GL_TEXTURE_2D, tex);
	unsigned char* image = SOIL_load_image("fire.jpeg", &width, &height, 0, SOIL_LOAD_RGBA);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);

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
	p.draw(tex);
	glPopMatrix();

}

/**
 * generateSmoke - Self Generates Smoke
 */
void generateSmoke()
{
	fluidSolver->resetInitialFields();

	int rowSize = fluidSolver->getRows();
	int colSize = fluidSolver->getCols();

	// Check Bounds
	if (xPos > 0 && xPos < rowSize - 1 && yPos > 0 && yPos < colSize - 1)
	{
		// Get the velocity values
		float xVel = 1.0f * x_velocity;
		float yVel = 1.0f * y_velocity;

		// Set the initial velocity
		fluidSolver->setInitVelocity(xPos, yPos, xVel, yVel);

		// Density Value
		float density = 1.0f;
		fluidSolver->setInitDensity(xPos, yPos, density);
	}

	// Set Velocities/Densities
	fluidSolver->addForce();
}

/**
 * getMouseInput
 */
void getMouseInput()
{
	// Reset initial velocities/densities
	fluidSolver->resetInitialFields();

	// Get Rows/Cols
	int rows = fluidSolver->getRows();
	int cols= fluidSolver->getCols();

	// Check if Left/Right Mouse Button Clicked
	if (_leftMouseButton || _rightMouseButton)
	{
		// Get Mouse position relative to grid
		int xPos = (int)((float)(_origMousePos[0]) / WIDTH * (rows));
		int yPos = (int)((float)(HEIGHT - _origMousePos[1]) / HEIGHT * (cols));

		// Bounds Check
		if (xPos > 0 && xPos < rows - 1 && yPos > 0 && yPos < cols - 1)
		{
			// Left Mouse Button
			if (_leftMouseButton)
			{
				// Get the velocity values
				float xVel = 1.0f * (_mousePos[0] - _origMousePos[0]);
				float yVel = 1.0f * (_origMousePos[1] - _mousePos[1]);

				// Set the initial velocity
				fluidSolver->setInitVelocity(xPos, yPos, xVel, yVel);
			}

			// Right Mouse Button
			if (_rightMouseButton)
			{
				// Density Value
				float density = 10.0f;
				fluidSolver->setInitDensity(xPos, yPos, density);
			}

			// Update original mouse position
			_origMousePos[0] = _mousePos[0];
			_origMousePos[1] = _mousePos[1];
		}

		// Set Velocities/Densities
		fluidSolver->addForce();
	}
}

/**
 * processKeys - Handles Callbacks for Keyboard Keys
 */
void processKeys(unsigned char key, int x, int y)
{
	switch (key)
	{
		// Change Display Mode
		case 'd':
		case 'D':
			displayMode = (displayMode + 1) % 2;
			break;

		case 'c':
		case 'C':
			fluidSolver->resetFields();
			x_velocity = 90;
			y_velocity = 120;
			break;
	}
}

/**
 * mouseDrag - Gets mouse drags information
 */
void mouseDrag(int x, int y)
{
	// Update mouse position
	_mousePos[0] = x;
	_mousePos[1] = y;
}


/**
 * mouseButton - Gets the State of the Mouse Button
 */
void mouseButton(int button, int state, int x, int y)
{
	// State of Mouse Button
	switch (button)
	{
		case GLUT_LEFT_BUTTON:
			_leftMouseButton = (state == GLUT_DOWN);
			break;
		case GLUT_MIDDLE_BUTTON:
			_middleMouseButton = (state == GLUT_DOWN);
			break;
		case GLUT_RIGHT_BUTTON:
			_rightMouseButton = (state == GLUT_DOWN);
			break;
	}

	// Update original mouse position
	_origMousePos[0] = x;
	_origMousePos[1] = y;

	// Update mouse position
	_mousePos[0] = x;
	_mousePos[1] = y;
}

/**
 * reshape
 */
void reshape(int w, int h)
{
	// Set Image Size
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// Set Camera Perspective
	glOrtho(0.0f, (float)(fluidSolver->getRows()) - BUFFER, 0.0f, (float)(fluidSolver->getCols()) - BUFFER, -LENGTH, LENGTH);

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
	//getMouseInput();
	generateSmoke();

	// Move the Velocity/Density forward 1 timestep
	fluidSolver->stepVelocity();
	fluidSolver->stepDensity();

	// Draw Density/velocity
	if (displayMode == 0) {
		drawDensity();
	}
	else {
		drawVelocity();
	}

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

	// Set Initial Position
	xPos = fluidSolver->getRows() / 2;
	yPos = fluidSolver->getCols() / 2;
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
	fluidSolver = new FluidSolver(100, 100, 1.0, 0.0f, 0.0f, 30);
	fluidSolver->resetFields();

	// Initialize GLUT
	glutInit(&argc, argv);

	// Request Color and Double Buffer
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

	// Set Window Size / Position
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(WIDTH, HEIGHT);

	// Create Window
	glutCreateWindow("Realistic Fire");

	

	// GLUT Callbacks
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(processKeys);
	glutMouseFunc(mouseButton);
	glutMotionFunc(mouseDrag);

	// Initialize States
	init();
	initTex();
	

	glutMainLoop();

	return 0;
}