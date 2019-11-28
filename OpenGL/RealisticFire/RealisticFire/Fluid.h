
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

	// Reset Functions
	void resetFields();
	void resetInitialFields();

	// Step 1: Add Forces
	void addForce();

	// Step 2: Advection
	void advection(float *value, float *value0, float *u, float *v);
	void advectVelocity(velocity* u1, velocity* u0);
	void advectDensity(float* d1, float* d0, velocity* u1);

	// Step 3: Diffusion
	void diffusion(float *value, float *value0, float rate);
	void diffuseVelocity(velocity* u1, velocity* u0, float rate);

	// Step 4: Projection
	void projection();

	// Index into Grid
	int idx(int i, int j);

	// Step Functions
	void stepVelocity();
	void stepDensity();

	void setBoundary(float *value, int flag);
	void setVelocityBoundary(velocity* u1);

	// Setters
	void setInitVelocity(int i, int j, float xVel, float yVel);
	void setInitDensity(int i, int j, float density);
	void setTimestep(float dt);
	void setViscosity(float visc);
	void setDiffusionRate(float diffRate);
	void setJacobiIterations(int jNum);

	// Getters
	point* getPoints();
	int getRows();
	int getCols();
	int getSize();
	velocity* getVelocities();
	float* getDensity();
	float getDens(int i, int j) { return (_d1[idx(i - 1, j - 1)] + _d1[idx(i, j - 1)] + _d1[idx(i - 1, j)] + _d1[idx(i, j)]) / 4.0f; }

	void swapVelocity(velocity* u1, velocity* u0);

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

	// Divergence
	float *divergence;

	// Pressure
	float *pressure;

	// Diffusion Properties
	float _viscosity;
	float _diffusionRate;

	// Grid Point Array
	point *p;
	
	// Velocity Fields
	velocity *_u1;
	velocity *_u0;

	// Density Fields
	float *_d1;
	float *_d0;

	float minX;
	float maxX;
	float minY;
	float maxY;
};

#endif // __FLUID_SOLVER_H__