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
void FluidSolver::init(int rowSize, int colSize)
{
	// Initialize Rows/Cols
	rows = rowSize;
	cols = colSize;
	size = rows * cols;

	// Initialize Velocities / Previous Velocity arrays
	u0 = new velocity[size];
	u1 = new velocity[size];

	// Initialize Grid array
	pos = new point[size];

	// Iterate over rows
	for (int i = 0; i < rows; ++i)
	{
		// Iterate over cols
		for (int j = 0; j < cols; ++j)
		{
			// Init grid points
			pos[idx(i, j)].x = (float)i + 0.5f;
			pos[idx(i, j)].y = (float)j + 0.5f;
		}
	}
}

/**
 * addSource
 */
void FluidSolver::addSource()
{
	// Iterate over rows
	for (int i = 1; i < rows-1; ++i)
	{
		// Iterate over cols
		for (int j = 1; j < cols-1; ++j)
		{
			// Get Index into Grid
			int c = idx(i, j);
			u1[c].x += u0[c].x;
			u1[c].y += u0[c].y;
		}
	}
}

/**
 * reset
 */
void FluidSolver::reset()
{
	// Iterate over Size of Grid
	for (int i = 0; i < size; i++)
	{
		// Initialize Velocity Values
		u1[i].x = 0.0;
		u1[i].y = 0.0;

		// Initialize Previous Velocity Values
		u0[i].x = 0.0;
		u0[i].y = 0.0;
	}
}

/**
 * resetInitialVelocities - Resets the initial velocities
 */
void FluidSolver::resetInitialVelocities()
{
	// Iterate over Grid
	for (int i = 0; i < size; i++)
	{
		u0[i].x = 0.0;
		u0[i].y = 0.0;
	}
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
//void FluidSolver::diffusion(velocity* u1, velocity* u0, float visc, float dt, int h)
//{
	// Iterate over rows
	//for (int i = 1; i < rows - 1; ++i)
	//{
		// Iterate over cols
		//for (int j = 1; j < cols - 1; ++j)
		//{
			// Get Cell Indices
			//int c = idx(i, j);
			//int c0 = idx(i-1, j); int c1 = idx(i+1, j);
			//int c2 = idx(i, j-1); int c3 = idx(i, j+1);

			// Get densities at the cells
			//velocity d0 = u1[c0];
			//velocity d1 = u1[c1];
			//velocity d2 = u1[c2];
			//velocity d3 = u1[c3];

			// Calculate Diffused values
			//u1[c].x = u0[c].x + (visc * dt) * (d0.x + d1.x + d2.x + d3.x - (4 * c)) / (h*h);
			//u1[c].y = u0[c].y + (visc * dt) * (d0.y + d1.y + d2.y + d3.y - (4 * c)) / (h*h);
		//}
	//}
//}


/**
 * setBoundary
 */
void FluidSolver::setBoundary()
{
}

/**
 * advection
 *
 * u1 - Velocity Field
 * u0 - Previous Velocity Field
 * dt - Time Step
 */
void FluidSolver::advection(velocity* u1, velocity* u0, float dt)
{
	// Iterate over Rows
	for (int i = 1; i < rows-1; ++i)
	{
		// Iterate over Cols
		for (int j = 1; j < cols-1; ++j)
		{
			// Get Index into Grid
			int c = idx(i, j);

			// Get Previous Position [xP = xG - dt * uG]
			point p0;
			p0.x = pos[c].x - dt * u1[c].x;
			p0.y = pos[c].y - dt * u1[c].y;

			// Neighbor Cell Indices
			int i0 = (int)(p0.x - 0.5); int i1 = i0 + 1;
			int j0 = (int)(p0.y - 0.5); int j1 = j0 + 1;

			// Neighbor Cells
			int c0 = idx(i0, j0); int c1 = idx(i1, j0);
			int c2 = idx(i0, j1); int c3 = idx(i1, j1);

			// Calculate new advected value via Bilinear Interpolation
			float b = pos[c2].y - p0.y;
			float t = 1.0 - b;
			float l = pos[c1].x - p0.x;
			float r = 1.0 - l;

			u1[c].x = t * (t*(l*u0[c0].x) + (r*u0[c1].x)) + b * ((l*u0[c2].x) + (r*u0[c3].x));
			u1[c].y = t * (t*(l*u0[c0].y) + (r*u0[c1].y)) + b * ((l*u0[c2].y) + (r*u0[c3].y));
		}
	}
}

/**
 * projection
 */
void FluidSolver::projection()
{

}

/**
 * idx - Gets the index based on i and j
 */
int FluidSolver::idx(int i, int j)
{
	return j * rows + i;
}

/**
 * stepVelocity - Moves the Velocity forward one time step
 *                The fields advected by itelf, the field diffuses
 *                due to viscous friction within the fluid. And
 *                finally the velocity is forced to conserve
 *                mass
 */
void FluidSolver::stepVelocity(float visc, float dt)
{
	// Swap previous velocities with current velocities
	SWAP(u0, u1);

	// Step 1: Advection
	advection(u1, u0, dt);

	// Step 2: Diffusion
	//diffusion(vels, vels0, visc, dt);

	// Step 3: Projection
	//projection();
}

/**
 * stepDiffusion
 */
void FluidSolver::stepDensity()
{

}

/**
 * Setters
 */
void FluidSolver::setInitVelocity(int i, int j, float xVel, float yVel)
{
	// Get Index
	int c = idx(i, j);

	// Update Velocity Value
	u0[c].x = xVel;
	u0[c].y = yVel;
}

/**
 * Getters
 */
int FluidSolver::getSize() { return size; }
int FluidSolver::getRows() { return rows; }
int FluidSolver::getCols() { return cols; }

point* FluidSolver::getPositions() { return pos; }
velocity* FluidSolver::getVelocities() { return u1; }