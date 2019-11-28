
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
	FluidSolver(int rows, int cols, float dt, float visc, float diffRate, int jNum);
	~FluidSolver();

	void init();

	// Reset Functions
	void resetFields();
	void resetInitialFields();

	void start() { running = 1; }
	void stop() { running = 0; }
	int isRunning() { return running; }

	// Step 1: Add Forces
	void addForce();


	void setBoundary(float *value, int flag);
	void projection();
	void advection(float *value, float *value0, float *u, float *v, int flag);

	// Index into Grid
	int idx(int i, int j);

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
	void setDiffusionRate(float diffRate);

	// Getters
	point* getPoints();
	int getRows();
	int getCols();
	int getSize();


	float getH() { return h; }
	float* getVX() { return vx; }
	float* getVY() { return vy; }
	float* getD() { return d1; }
	float getDens(int i, int j) { return (d1[idx(i - 1, j - 1)] + d1[idx(i, j - 1)] + d1[idx(i - 1, j)] + d1[idx(i, j)]) / 4.0f; }

private:

	// Total Grid Size
	int _size;

	// X/Y Dimensions
	int _rows;
	int _cols;

	// TimeStep
	float _dt;

	// Number of Jacobi Iterations
	int _jNum;

	// Diffusion Properties
	float _viscosity;
	float _diffusionRate;

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
	float minX;
	float maxX;
	float minY;
	float maxY;

	//params
	int running;

	float *div;
	float *p;
};

#endif // __FLUID_SOLVER_H__