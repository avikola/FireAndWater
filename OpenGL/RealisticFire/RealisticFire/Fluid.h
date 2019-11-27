
#ifndef __FLUID_SOLVER_H__
#define __FLUID_SOLVER_H__

// Point Object
struct point
{
	float x;
	float y;
};

// Velocity Object
struct velocity
{
	float x;
	float y;
};

class FluidSolver
{
public:

	// Constructor/Destructor
	FluidSolver(int rows, int cols, float dt, float visc);
	~FluidSolver();

	void init();

	// Reset Functions
	void resetFields();
	void resetInitialFields();

	void start() { running = 1; }
	void stop() { running = 0; }
	int isRunning() { return running; }

	//animation
	void setBoundary(float *value, int flag);
	void projection();
	void advection(float *value, float *value0, float *u, float *v, int flag);
	void vortConfinement();
	void addSource();

	// Diffusion Functions
	void diffusion(float *value, float *value0, float rate, int flag);
	void diffuseVelocity(velocity *u1, velocity *u0, float rate, int flag);

	// Step Functions
	void stepVelocity();
	void stepDensity();

	// Setters
	void setInitVelocity(int i, int j, float xVel, float yVel);
	void setInitDensity(int i, int j, float density);
	void setTimestep(float dt);
	void setViscosity(float visc);

	// Getters
	point* getPoints();
	int getRows();
	int getColSize() { return _cols; }
	int getTotSize() { return _size; }
	float getH() { return h; }
	float getSimSizeX() { return simSizeX; }
	float getSimSizeY() { return simSizeY; }
	float* getVX() { return vx; }
	float* getVY() { return vy; }
	float* getD() { return d1; }
	float getDens(int i, int j) { return (d1[cIdx(i - 1, j - 1)] + d1[cIdx(i, j - 1)] + d1[cIdx(i - 1, j)] + d1[cIdx(i, j)]) / 4.0f; }

	int cIdx(int i, int j) { return j * _rows + i; }

private:

	// Total Grid Size
	int _size;

	// X/Y Dimensions
	int _rows;
	int _cols;

	// TimeStep
	float _dt;

	// Viscosity
	float _viscosity;

	// Point Array
	point *pt;
	
	// Velocity Fields
	float *vx;
	float *vy;
	float *v0x;
	float *v0y;
	velocity *_u1;
	velocity *_u0;

	// Density Field
	float *d1;
	float *d0;

	float h;
	float simSizeX;
	float simSizeY;
	float minX;
	float maxX;
	float minY;
	float maxY;

	//params
	int running;

	float diff;
	float vorticity;
	float *div;
	float *p;
	//vorticity confinement
	float *vort;
	float *absVort;
	float *gradVortX;
	float *gradVortY;
	float *lenGrad;
	float *vcfx;
	float *vcfy;
};

#endif // __FLUID_SOLVER_H__