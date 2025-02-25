// Standard Libraries
#include <time.h>
#include <string>

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

// Initial Velocities + Direction
int radius;
int angle;
float radians;
float x_velocity;	// x velocity
float y_velocity;	// y velocity

// Velocity Scalers
float x_scaler;
float y_scaler;

// Cursor Velocity Scalers
float x_cursor_velocity;
float y_cursor_velocity;

// X Pos, Y Pos
static int xPos;
static int yPos;

// Grid resolution
int x_grid_width = 120;
int y_grid_height = 120;

// Smoke Viscosity
float density_factor;

// Smoke Amount
float density = 1.0f;

// Hide Text
bool text_visibility = true;

// Determine Mode - Smoke / Smoke + Fire
bool SMOKE_ONLY = true;

/**
 * drawVelocity
 */
void drawVelocity()
{	
	// Get Points
	point* p = fluidSolver->getPoints();

	// Get Velocities
	velocity* v = fluidSolver->getVelocities();

	// Line Color / Style
	glColor3f(0.0f, 0.0f, 1.0f);
	glLineWidth(1.0f);

	// Draw Velocity Lines
	glBegin(GL_LINES);
		for (int i = 0; i < fluidSolver->getSize(); ++i)
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
	// Get Rows / Cols
	int rowSize = fluidSolver->getRows();
	int colSize = fluidSolver->getCols();
	
	// Draw Fluid Density
	glBegin(GL_QUADS);
		for (int i = 1; i <= rowSize - 2; ++i)
		{
			float x = (float)i;
			for (int j = 1; j <= colSize - 2; ++j)
			{
				float y = (float)j;
				
				float d00 = fluidSolver->getDensity(i, j, density_factor);
				float d01 = fluidSolver->getDensity(i, j + 1, density_factor);
				float d10 = fluidSolver->getDensity(i + 1, j, density_factor);
				float d11 = fluidSolver->getDensity(i + 1, j + 1, density_factor);


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
void initTex() 
{
	// Image dimensions
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

	// Draw particles
	glPushMatrix();

	p.advance(DELTA);
	p.delete_particle();
	p.draw(tex);

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
	int rows = fluidSolver->getRows();
	int cols = fluidSolver->getCols();

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
				float xVel = x_cursor_velocity * (_mousePos[0] - _origMousePos[0]);
				float yVel = y_cursor_velocity * (_origMousePos[1] - _mousePos[1]);

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
 * generateSmoke - Self Generates Smoke
 */
void generateSmoke()
{
	// Get input from Mouse
	getMouseInput();

	// Generate Smoke

	fluidSolver->resetInitialFields();

	int rowSize = fluidSolver->getRows();
	int colSize = fluidSolver->getCols();

	// Check Bounds
	if (xPos > 0 && xPos < rowSize - 1 && yPos > 0 && yPos < colSize - 1)
	{
		// Get the velocity values
		float xVel = x_scaler * x_velocity;
		float yVel = y_scaler * y_velocity;

		// Set the initial velocity
		fluidSolver->setInitVelocity(xPos, yPos, xVel, yVel);
		fluidSolver->setInitDensity(xPos, yPos, density);
	}

	// Set Velocities/Densities
	fluidSolver->addForce();

	// Move the Velocity/Density forward 1 timestep
	fluidSolver->stepVelocity();
	fluidSolver->stepDensity();

	// Draw Density/velocity
	if (displayMode == 0)
		drawDensity();
	else
	{
		drawDensity();
		drawVelocity();
	}

}

// Change Smoke Position.
void smokeReposition(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_UP:
		if (yPos < y_grid_height - 5)
			yPos += 5;
		break;
	case GLUT_KEY_DOWN:
		if (yPos > 5)
			yPos -= 5;
		break;
	case GLUT_KEY_LEFT:
		if (xPos > 5)
			xPos -= 5;
		break;
	case GLUT_KEY_RIGHT:
		if (xPos < x_grid_width - 5)
			xPos += 5;
		break;
	}
	glutPostRedisplay();
}

void setVelDirection()
{	
	radians = angle * 0.0174532925f;
	x_velocity = radius * cos(radians);
	y_velocity = radius * sin(radians);
}

// Set various default parametric values
void setValues()
{
	x_scaler = 1.0;
	y_scaler = 1.0;
	x_cursor_velocity = 3.0;
	y_cursor_velocity = 3.0;
	density_factor = 4.5;
	radius = 40;
	angle = 90;
	setVelDirection();
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
		case 'c':	// Clear screen
		case 'C':
			fluidSolver->resetFields();
			break;
		case 'r':	// Reset Smoke
		case 'R':
			fluidSolver->resetFields();
			
			yPos = (fluidSolver->getCols() / 4) - 10;
			(SMOKE_ONLY) ? xPos = fluidSolver->getRows() / 2 : xPos = (fluidSolver->getRows() / 2) + 10;
			
			setValues();

			break;
		case 'h':	// Hide text overlay
		case 'H':
			text_visibility = !text_visibility;
			break;
		case 'k':	// Increase velocity
		case 'K':
			y_scaler += 0.5;
			break;
		case 'j':	// Decrease velocity
		case'J':
			y_scaler -= 0.5;
			break;
		case 'p':	// Rotate right
		case 'P':
			if (angle == 0)
				angle = 360;
			angle -= 5;

			setVelDirection();
			break;
		case 'o':	// Rotate left
		case 'O':
			if (angle == 360)
				angle = 0;
			angle += 5;

			setVelDirection();
			break;
		case 'v':	// Decrease cursor velocity
		case 'V':
			x_cursor_velocity -= 0.5;
			y_cursor_velocity -= 0.5;
			break;
		case 'b':	// Increase cursor velocity
		case 'B':
			x_cursor_velocity += 0.5;
			y_cursor_velocity += 0.5;
			break;
		case 'f':	// Decrease density factor
		case 'F':
			if (density_factor != 1.5)
				density_factor -= 0.5;
			break;
		case 'g':	// Increase density factor
		case 'G':
			density_factor += 0.5;
			break;
		case '1':	// Smoke Only mode
			SMOKE_ONLY = true;
			text_visibility = true;

			xPos = fluidSolver->getRows() / 2;
			break;
		case '2':	// Smoke + Fire mode
			SMOKE_ONLY = false;

			xPos = (fluidSolver->getRows() / 2) + 10;
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
	glutReshapeWindow(WIDTH, HEIGHT);

	// Set Image Size
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// Set Camera Perspective
	glOrtho(BUFFER, (float)(fluidSolver->getRows()) - BUFFER, BUFFER, (float)(fluidSolver->getCols()) - BUFFER, -LENGTH, LENGTH);

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
	// Update screen
	glutPostRedisplay();
}

// Function to display text.
void drawStrings(const char* str, int len, int x, int y)
{
	glMatrixMode(GL_PROJECTION);
	double* mat = new double[16];
	glGetDoublev(GL_PROJECTION_MATRIX, mat);
	glLoadIdentity();
	glOrtho(0, WIDTH, 0, HEIGHT, -5, 5);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glPushMatrix();
	glLoadIdentity();
	glColor3f(1.0f, 0.2f , 0.0f);
	glRasterPos2i(x, y);
	for (int i = 0; i < len; ++i)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, (int)str[i]);
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glLoadMatrixd(mat);
	glMatrixMode(GL_MODELVIEW);
}

// Text Function
void setText()
{
	// Draw Text
	if (text_visibility)
	{
		string intro = "Smoke Controls:";
		drawStrings(intro.data(), intro.size(), 5, HEIGHT - 20);

		string arrow_key_hint = "Arrow Keys to Change Position";
		drawStrings(arrow_key_hint.data(), arrow_key_hint.size(), 5, HEIGHT - 40);

		string current_pos = "Position: (" + to_string(xPos) + ", " + to_string(yPos) + ")";
		drawStrings(current_pos.data(), current_pos.size(), 5, HEIGHT - 60);

		string angles = "'O' and 'P' to adjust rotation";
		drawStrings(angles.data(), angles.size(), 5, HEIGHT - 100);

		float temp = ceilf(y_scaler * 100) / 100;
		string v_factor = "Velocity Scale: " + to_string(temp);
		drawStrings(v_factor.data(), v_factor.size(), 5, HEIGHT - 140);
		v_factor = "'J' and 'K' to adjust";
		drawStrings(v_factor.data(), v_factor.size(), 5, HEIGHT - 160);

		string curse = "Cursor Velocity Scale: " + to_string(x_cursor_velocity);
		drawStrings(curse.data(), curse.size(), 5, HEIGHT - 200);
		curse = "'V' and 'B' to adjust";
		drawStrings(curse.data(), curse.size(), 5, HEIGHT - 220);
		curse = "Right-click to create smoke";
		drawStrings(curse.data(), curse.size(), 5, HEIGHT - 260);

		string dense = "Density: " + to_string(density_factor);
		drawStrings(dense.data(), dense.size(), 5, HEIGHT - 300);
		dense = "'F' and 'G' to adjust";
		drawStrings(dense.data(), dense.size(), 5, HEIGHT - 320);

		string hider = "'C' to clear screen";
		drawStrings(hider.data(), hider.size(), 5, HEIGHT - 360);
		hider = "'R' to reset";
		drawStrings(hider.data(), hider.size(), 5, HEIGHT - 380);
		hider = "'H' to hide text";
		drawStrings(hider.data(), hider.size(), 5, HEIGHT - 400);
	}
}

/**
 * display
 */
void display()
{
	if (SMOKE_ONLY)
		generateSmoke();

	else
	{
		generateSmoke();
		drawFire();
	}

	setText();

	// Double Buffer Flush
	glutSwapBuffers();

	glutPostRedisplay();
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

	// Enable transparency
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Set background color [Black]
	glClearColor(0.0, 0.0, 0.0, 0.0); 

	// Set Initial Position
	xPos = fluidSolver->getRows() / 2;
	yPos = (fluidSolver->getCols() / 4) - 10 ;

	// Default Values
	setValues();
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
	fluidSolver = new FluidSolver(x_grid_width, y_grid_height, 0.4f, 0.0f, 0.0f, 30);
	fluidSolver->resetFields();

	// Initialize GLUT
	glutInit(&argc, argv);

	// Request Color and Double Buffer
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

	// Set Window Size / Position
	glutInitWindowPosition(150, 10);
	glutInitWindowSize(WIDTH, HEIGHT);

	// Create Window
	glutCreateWindow("Realistic Fire + Smoke");

	// GLUT Callbacks
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(processKeys);
	glutMouseFunc(mouseButton);
	glutMotionFunc(mouseDrag);
	glutSpecialFunc(smokeReposition);

	// Initialize States
	init();
	initTex();
	
	glutMainLoop();

	return 0;
}