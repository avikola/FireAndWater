#include "FluidSolver.h"

/**
 * Constructor
 */
FluidSolver::FluidSolver()
{
}

/**
 * Destructor
 */
FluidSolver::~FluidSolver()
{
}

/**
 * init
 */
void FluidSolver::init(int rows, int cols)
{
	// Initialize Rows/Cols
	rows = rows; cols = cols;

	size = rows * cols;

	// Initialize Timestep
	t = 1.0;
}

/**
 * diffusion - Calculates the effects of diffusion on the fluid
 *
 * d0 - Density grid before diffusion
 * d  - Density grid after diffusion
 * dt - Time Step
 * r  - Diffusion Rate
 * h  - Grid Spacing
 */
void FluidSolver::diffusion(float *val, float *val0, float dt, float r, float h)
{
	// Iterate over rows
	for (int i = 1; i < rows - 1; ++i)
	{
		// Iterate over cols
		for (int j = 1; j < cols - 1; ++j)
		{
			// Get Cell Indices
			int c = cellIndex(i, j);
			int c0 = cellIndex(i-1, j);
			int c1 = cellIndex(i+1, j);
			int c2 = cellIndex(i, j-1);
			int c3 = cellIndex(i, j+1);

			// Get densities at the cells
			int d0 = val[c0];
			int d1 = val[c1];
			int d2 = val[c2];
			int d3 = val[c3];

			// Calculate Diffused Density
			val[c] = val0[c] + r * dt * (d0 + d1 + d2 + d3 - (4 * c)) / (h*h);
		}
	}
}


/**
 * setBoundary
 */
void FluidSolver::setBoundary()
{
	for (int i = 1; i < rows - 1; i++)
	{

	}
}

/**
 * advection
 *
 * u - Velocity Field (X Axis)
 * v - Velocity Field (Y Axis)
 * dt - Time Step
 */
void FluidSolver::advection(float *val0, float *val, float *u, float *v, float dt)
{
	// Iterate over cell rows
	for (int i=1; i<rows-1; i++)
	{
		// Iterate over cell cols
		for (int j=1; j<cols-1; j++)
		{
			// Get Cell
			int c = cellIndex(i, j);

			// Get Previous Position
			float x0 = px[c] - u[c] * dt;
			float y0 = py[c] - v[c] * dt;

			// Neighbor Cell Indices
			int i0 = (int)x0; int i1 = i0+1;
			int j0 = (int)y0; int j1 = j0+1;

			// Neighbor Cells
			int c0 = cellIndex(i0, j0); int c1 = cellIndex(i1, j0);
			int c2 = cellIndex(i0, j1); int c3 = cellIndex(i1, j1);

			float t;
			float l;
			float r;
			float b;

			// Calculate new advected value via Bilinear Interpolation
			val[c] = t * (t*(l*val0[c0]) + (r*val0[c1])) +  b*((l*val0[c2]) + (r*val0[c3]));
		}
	}

	// Update Boundarys of Fluid Container
	setBoundary();
}

void FluidSolver::projection()
{

}

/**
 * cellIndex - Gets the Cell Index based on i and j
 */
int FluidSolver::cellIndex(int i, int j)
{
	return j * rows + i;
}

/**
 * stepVelocity
 */
void FluidSolver::stepVelocity()
{
	// Step 1: Advection
	advection(vX, v0X, v0X, v0Y, t);
	advection(vY, v0Y, v0X, v0Y, t);

	// Step 2: Diffusion
	diffusion(vX, v0X, t, 0, 1);
	diffusion(vY, v0Y, t, 0, 1);

	// Step 3: Projection
	projection();
}

/**
 * stepDiffusion
 */
void FluidSolver::stepDensity()
{

}

/**
 * Getters
 */
float* FluidSolver::getVelocityX() { return vX; }
float* FluidSolver::getVelocityY() { return vY; }
point* FluidSolver::getPositions() { return pos; }
