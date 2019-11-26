// System Libraries
#include <time.h>
#include <iostream>
#include <stdlib.h>
#include <windows.h>
#include <string>
#include <experimental/filesystem>

// Headers
#include "core.h"
#include "controller.h"
#include "fluid.h"
#include "particle_system.h"
#include "RealisticFire.h"
#include "constants.h"
#include "cloud.h"
#include "imageIO.h"
//#include "SOIL.h"

// Namespaces
using namespace std;

#define WINDOW_NAME		 "Smoke3D"

Controller* g_controller = NULL;
Fluid *fluid;

particle_system p(NUMBER_OF_PARTICLES);

// Window Width/Height
float winWidth = 400, winHeight = 400;

// Current Mouse Position
int _mousePos[2] = { 0, 0 };

// Enabled Indicators
bool _translateEnabled;

//  Mouse Buttons
int _leftMouseButton = 0;
int _middleMouseButton = 0;
int _rightMouseButton = 0;

// Number of Textures
GLuint groundTexHandle;

// State of the World Transformations
float _theta[3] = { 0, 0, 0 };
float _translate[3] = { 0, 0, 0 };
float _scale[3] = { 1, 1, 1 };

// Initialize Display List Index
GLuint _listIndex;

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

		// Translate Case
		case 't':
		case 'T':
			// Check if Translate already Enabled
			if (_translateEnabled)
			{
				// Disable Translate Mode
				_translateEnabled = false;
			}
			else
			{
				// Enable Translate Mode
				_translateEnabled = true;
			}

			// Update Mouse Position
			_mousePos[0] = x;
			_mousePos[1] = y;
			break;


		// Z Key (Moves Camera in -Y Direction)
		case 'z':
		case 'Z':
			_translate[1] -= 1.0;
			break;

		// X Key (Moves Camera in Y Direction)
		case 'x':
		case 'X':
			_translate[1] += 1.0;
			break;

		case 27:
			exit(0);
			break;
	}
	glutPostRedisplay();
}

/**
 * processSpecialKeys - Handles Callbacks for Arrow Keys
 */
void processSpecialKeys(int key, int x, int y)
{
	switch (key)
	{
		// Up Key (Moves Camera in -X Direction)
	case GLUT_KEY_UP:
		_translate[0] -= 1.0;
		break;

		// Down Key (Moves Camera in X Direction)
	case GLUT_KEY_DOWN:
		_translate[0] += 1.0;
		break;

		// Left Key (Moves Camera in +Z Direction)
	case GLUT_KEY_LEFT:
		_translate[2] += 1.0;
		break;

		// Right Key (Moves Camera in -Z Direction)
	case GLUT_KEY_RIGHT:
		_translate[2] -= 1.0;
		break;
	}

	// Make the screen update
	glutPostRedisplay();
}

/**
 * drawGround - Initializes the Ground Texture. It then
 *              stores the Vertices into a Display List
 *              be later displayed
 */
void drawGround()
{
	// Compile the Ground into the Display List
	glNewList(_listIndex, GL_COMPILE);

	// Activate Ground Texture
	glEnable(GL_TEXTURE_2D);

	// Load Ground Texture
	glBindTexture(GL_TEXTURE_2D, groundTexHandle);

	// Initialize Texture Properties
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // Repeat Pattern in S
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // Repeat Pattern in T
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); // Bilinear Interpolation
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); // Bilinear Interpolation

   // Draw Ground
	glBegin(GL_QUADS);
	glTexCoord2f(0.0, 1.0f); glVertex3f(200.0f, -10.0, 200.0f);     // Vertex 1
	glTexCoord2f(0.0f, 0.0f); glVertex3f(200.0f, -10.0, -200.0f);   // Vertex 2
	glTexCoord2f(1.0f, 1.0f); glVertex3f(-200.0f, -10.0, -200.0f);  // Vertex 3
	glTexCoord2f(1.0f, 0.0f); glVertex3f(-200.0f, -10.0, 200.0f);   // Vertex 4
	glEnd();

	// Disable Ground Texture
	glDisable(GL_TEXTURE_2D);

	// End Compilation of the Ground
	glEndList();
}

int initTexture(const char * imageFilename, GLuint textureHandle)
{
	// read the texture image
	ImageIO img;
	ImageIO::fileFormatType imgFormat;
	ImageIO::errorType err = img.load(imageFilename, &imgFormat);

	if (err != ImageIO::OK)
	{
		printf("Loading texture from %s failed.\n", imageFilename);
		return -1;
	}

	// check that the number of bytes is a multiple of 4
	if (img.getWidth() * img.getBytesPerPixel() % 4)
	{
		printf("Error (%s): The width*numChannels in the loaded image must be a multiple of 4.\n", imageFilename);
		return -1;
	}

	// allocate space for an array of pixels
	int width = img.getWidth();
	int height = img.getHeight();
	unsigned char * pixelsRGBA = new unsigned char[4 * width * height]; // we will use 4 bytes per pixel, i.e., RGBA

	// fill the pixelsRGBA array with the image pixels
	memset(pixelsRGBA, 0, 4 * width * height); // set all bytes to 0
	for (int h = 0; h < height; h++)
		for (int w = 0; w < width; w++)
		{
			// assign some default byte values (for the case where img.getBytesPerPixel() < 4)
			pixelsRGBA[4 * (h * width + w) + 0] = 0; // red
			pixelsRGBA[4 * (h * width + w) + 1] = 0; // green
			pixelsRGBA[4 * (h * width + w) + 2] = 0; // blue
			pixelsRGBA[4 * (h * width + w) + 3] = 255; // alpha channel; fully opaque

			// set the RGBA channels, based on the loaded image
			int numChannels = img.getBytesPerPixel();
			for (int c = 0; c < numChannels; c++) // only set as many channels as are available in the loaded image; the rest get the default value
				pixelsRGBA[4 * (h * width + w) + c] = img.getPixel(w, h, c);
		}

	// bind the texture
	glBindTexture(GL_TEXTURE_2D, textureHandle);

	// initialize the texture
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixelsRGBA);

	// generate the mipmaps for this texture
	//glGenerateMipmap(GL_TEXTURE_2D);

	// set the texture parameters
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	// query support for anisotropic texture filtering
	GLfloat fLargest;
	glGetFloatv(GL_MAX_TEXTURE_MAX_ANISOTROPY_EXT, &fLargest);
	printf("Max available anisotropic samples: %f\n", fLargest);
	// set anisotropic texture filtering
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 0.5f * fLargest);

	// query for any errors
	GLenum errCode = glGetError();
	if (errCode != 0)
	{
		printf("Texture initialization error. Error code: %d.\n", errCode);
		return -1;
	}

	// de-allocate the pixel array -- it is no longer needed
	delete[] pixelsRGBA;

	return 0;
}

//Initializes 3D rendering
void init()
{
	// Make big points and wide lines
	glPointSize(1);

	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

	// Initialize Keyboard Indicators
	_translateEnabled = false;

	//Enable transparency
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Initialize Display List
	_listIndex = glGenLists(1);

	// Initialize Textures
	glGenTextures(1, &groundTexHandle);
	int code = initTexture("ground.jpg", groundTexHandle);
	if (code != 0)
	{
		printf("Error loading the ground texture image.\n");
		exit(EXIT_FAILURE);
	}

	// Initialize Objects
	drawGround();
}

/*
 * reshape - Called every time window is resized
 *           to update the projection matrix, and
 *           to preserve aspect ratio
 */
void reshape(int w, int h)
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

	// Initialize Aspect
	GLfloat aspect = (GLfloat)w / (GLfloat)h;

	// Setup image size
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// Setup camera
	gluPerspective(60.0, aspect, 0.1, 100.0);

	glViewport(0, 0, winWidth, winHeight);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//glOrtho(-LENGTH, LENGTH, -LENGTH, LENGTH, -LENGTH, LENGTH);

	glMatrixMode(GL_MODELVIEW);
}


//
// drawFire - Draws the Fire
//
void drawFire()
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
		glVertex3f(50, -5, 50);
		glVertex3f(50, -5, -50);
		glVertex3f(-50, -5, -50);
		glVertex3f(-50, -5, 50);
	glEnd();

	glutSwapBuffers();
	glutPostRedisplay();
}

/**
 * display
 */
void display()
{
	// Clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Reset transformation
	glLoadIdentity();

	// Use Texture Color Directly
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	// Transformations
	glScalef(_scale[0], _scale[1], _scale[2]);
	glRotatef(_theta[0], 1, 0, 0);
	glRotatef(_theta[1], 0, 1, 0);
	glRotatef(_theta[2], 0, 0, 1);
	glTranslatef(_translate[0], _translate[1], _translate[2]);

	// Draw Functions
	glCallList(_listIndex);   // Ground
	drawFire();               // Fire

	//fluid->SimulateStep();
	//fluid->Show();
}

/**
 * mousebutton
 */
void mousebutton(int button, int state, int x, int y)
{
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

	switch (glutGetModifiers())
	{
	case GLUT_ACTIVE_CTRL:
		g_ControlState = TRANSLATE;
		break;

	case GLUT_ACTIVE_SHIFT:
		g_ControlState = SCALE;
		break;

	default:
		if (_translateEnabled)
		{
			g_ControlState = TRANSLATE;
		}
		else
		{
			g_ControlState = ROTATE;
		}
		break;
	}

	_mousePos[0] = x;
	_mousePos[1] = y;
}

/**
 * mouseDrag - Converts mouse drags into information
 *             about rotation/translation/scaling
 */
void mousedrag(int x, int y)
{
	// Update Change in Mouse Position
	int vMouseDelta[2] = { x - _mousePos[0], y - _mousePos[1] };

	// Check Control State
	switch (g_ControlState)
	{
		// Translate Transformation
	case TRANSLATE:
		if (_leftMouseButton)
		{
			_translate[0] += 0.01 * vMouseDelta[0];
			_translate[1] -= 0.01 * vMouseDelta[1];
		}
		else if (_middleMouseButton)
		{
			_translate[2] += 0.01 * vMouseDelta[1];
		}
		break;

		// Rotate Transformation
	case ROTATE:
		if (_leftMouseButton)
		{
			_theta[0] += vMouseDelta[1];
			_theta[1] += vMouseDelta[0];
		}
		else if (_middleMouseButton)
		{
			_theta[2] += vMouseDelta[1];
		}
		break;

		// Scale Transformation
	case SCALE:
		if (_leftMouseButton)
		{
			_scale[0] *= 1.0 + vMouseDelta[0] * 0.01;
			_scale[1] *= 1.0 - vMouseDelta[1] * 0.01;
		}
		else if (_middleMouseButton)
		{
			_scale[2] *= 1.0 - vMouseDelta[1] * 0.01;
		}
		break;
	}

	// Update Current Mouse Position
	_mousePos[0] = x;
	_mousePos[1] = y;
}


//Handle mouse movement
void mouse_movement(int x, int y)
{
	float ww_ratio = float(x) / winWidth;
	float wh_ratio = float(y) / winHeight;

	p.set_gravity(vec3d((2 * ww_ratio - 1)*LENGTH, (1 - 2 * wh_ratio)*LENGTH, 0));
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

//
// main
//
int main(int argc, char** argv)
{
	//srand(time(0));
	//p.set_gravity();
	
	//::g_controller = new Controller(argc, argv, WINDOW_NAME);


	//::g_controller->InitCamera();
	//::g_controller->RegisterObject(object);
	//::g_controller->BeginLoop();

	//Initialize GLUT
	glutInit(&argc, argv);

	// Initialize Fluid
	fluid = new Fluid;

	// Request double buffer, depth, and color
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

	// Set Window Size and Position
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(WIN_WIDTH, WIN_HEIGHT);

	//Create the window and initialize OpenGL
	glutCreateWindow("Realistic Fire");

	// GLUT Callbacks
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutKeyboardFunc(handle_keypress);
	glutReshapeFunc(reshape);
	glutMouseFunc(mousebutton);
	glutMotionFunc(mousedrag);
	glutPassiveMotionFunc(mouse_movement);
	glutSpecialFunc(processSpecialKeys);

	// Initialize States
	init();

	glutMainLoop();

	return 0;
}