
#ifndef __FLUID_SOLVER_H__
#define __FLUID_SOLVER_H__

// Point Object
struct point
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

	// Step 3: Diffusion
	void diffusion(float *value, float *value0, float rate);

	// Step 4: Projection
	void projection();

	// Index into Grid
	int idx(int i, int j);

	// Step Functions
	void stepVelocity();
	void stepDensity();

	void setBoundary(float *value, int flag);

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

	float* getVX() { return ux; }
	float* getVY() { return uy; }
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
	float *ux;
	float *uy;
	float *v0x;
	float *v0y;

	// Density Field
	float *d1;
	float *d0;

	float minX;
	float maxX;
	float minY;
	float maxY;
};

#endif // __FLUID_SOLVER_H__