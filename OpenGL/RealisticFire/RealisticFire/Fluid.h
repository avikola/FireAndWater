
#ifndef __FLUID_SOLVER_H__
#define __FLUID_SOLVER_H__

#define SWAP(value0,value) {float *tmp=value0;value0=value;value=tmp;}

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
	void advectVelocity(velocity* u1, velocity* u0);
	void advectDensity(float* d1, float* d0, velocity* u1);

	// Step 3: Diffusion
	void diffusion(float *value, float *value0, float rate);
	void diffuseVelocity(velocity* u1, velocity* u0, float rate);

	// Step 4: Projection
	void projection();
	void computePressure();

	// Index into Grid
	int idx(int i, int j);

	// Step Functions
	void stepVelocity();
	void stepDensity();

	// Swap Functions
	void swapVelocity(velocity* u1, velocity* u0);

	// Boundary Cases
	void setPressureBoundary(float *pressure);
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
	float* getDensities();
	float getDensity(int i, int j);

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
	float *_pressure;

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