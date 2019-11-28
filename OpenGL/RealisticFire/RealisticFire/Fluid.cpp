// Standard Libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Headers
#include "Fluid.h"

#define SWAP(value0,value) {float *tmp=value0;value0=value;value=tmp;}

/**
 * Constructor
 *
 * rows     - Rows in Grid
 * cols     - Cols in Grid
 * dt       - TimeStep
 * visc     - Viscosity
 * diffRate - Diffusion Rate
 * jNum     - Jacobian Iterations
 */
FluidSolver::FluidSolver(int rows, int cols, float dt, float visc, float diffRate, int jNum)
{
	// Set Dimension Sizes
	_rows = rows;
	_cols = cols;

	// Set Grid Size
	_size = _rows * _cols;

	// Initialize TimeStep
	_dt = dt;

	// Initialize Number of Jacobi Iterations
	_jNum = jNum;

	// Initialize Diffusion properties
	_viscosity = visc;
	_diffusionRate = diffRate;

	// Initialize Point Array
	p = new point[_size];

	// Initialize Velocity Arrays
	ux = new float[_size];
	uy = new float[_size];
	v0x = new float[_size];
	v0y = new float[_size];

	// Initialize Density Arrays
	d1 = new float[_size];
	d0 = new float[_size];

	// Initialize Divergence
	divergence = new float[_size];

	// Initialize Pressure
	pressure = new float[_size];

	minX = 1.0f;
	maxX = _rows - 1.0f;
	minY = 1.0f;
	maxY = _cols - 1.0f;

	for (int i = 0; i < _rows; i++)
	{
		for (int j = 0; j < _cols; j++)
		{
			p[idx(i, j)].x = (float)i + 0.5f;
			p[idx(i, j)].y = (float)j + 0.5f;
		}
	}
}

/**
 * Destructor
 */
FluidSolver::~FluidSolver(){}

/**
 * addForce - Adds the Initial Force to velocity/density
 */
void FluidSolver::addForce()
{
	// Iterate over the Grid
	for (int i = 1; i <= _size; i++)
	{
		// Set Velocity
		ux[i] += v0x[i];
		uy[i] += v0y[i];

		// Set Density
		d1[i] += d0[i];
	}
}

/**
 * advection
 */
void FluidSolver::advection(float *value, float *value0, float *u, float *v)
{
	for (int i = 0; i < _rows; i++)
	{
		for (int j = 0; j < _cols; j++)
		{
			// Get grid cell
			int c = idx(i, j);

			float oldX = p[c].x - u[c] * _dt;
			float oldY = p[c].y - v[c] * _dt;

			if (oldX < minX) oldX = minX;
			if (oldX > maxX) oldX = maxX;
			if (oldY < minY) oldY = minY;
			if (oldY > maxY) oldY = maxY;

			int i0 = (int)(oldX - 0.5f);
			int j0 = (int)(oldY - 0.5f);
			int i1 = (i0 + 1) % _rows;
			int j1 = (j0 + 1) % _cols;

			// Perform Bilinear Interpolation
			float wL = p[idx(i1, j0)].x - oldX;
			float wR = 1.0f - wL;
			float wB = p[idx(i0, j1)].y - oldY;
			float wT = 1.0f - wB;

			value[idx(i, j)] = wB * (wL*value0[idx(i0, j0)] + wR * value0[idx(i1, j0)]) +
							   wT * (wL*value0[idx(i0, j1)] + wR * value0[idx(i1, j1)]);
		}
	}
}

/**
 * diffusion
 */
void FluidSolver::diffusion(float *v, float *v0, float rate)
{
	for (int i = 0; i < _size; i++)
	{
		v[i] = 0.0f;
	}
	float a = rate * _dt;

	// Jacobi Method Iteration
	for (int k = 0; k < _jNum; k++)
	{
		// Iterate over rows in grid
		for (int i = 0; i < _rows; i++)
		{
			// Iterate over cols in grid
			for (int j = 0; j < _cols; j++)
			{
				// Get grid cell
				int c = idx(i, j);

				// Check for Boundary cases
				int i0 = (_rows + ((i-1)%_rows)) + _rows;
				int j0 = (_cols + ((j-1)%_cols)) + _cols;
				int i1 = (i + 1) % _rows;
				int j1 = (j + 1) % _cols;
				
				// Neighbor grid Cells
				int c0 = idx(i0, j); int c1 = idx(i1, j);
				int c2 = idx(i, j0); int c3 = idx(i, j1);

				v[c] = (v0[c] + a * (v[c1] + v[c0] + v[c3] + v[c2])) / (4.0f*a + 1.0f);
			}
		}
	}
}

/**
 * projection
 */
void FluidSolver::projection()
{
	for (int i = 0; i < _rows; i++)
	{
		for (int j = 0; j < _cols; j++)
		{
			// Get grid cell
			int c = idx(i, j);

			// Neighbor Cells
			int c0 = idx(i - 1, j); int c1 = idx(i + 1, j);
			int c2 = idx(i, j - 1); int c3 = idx(i, j + 1);

			// Calculate Divergence
			divergence[c] = 0.5f * (ux[c1] - ux[c0] + uy[c3] - uy[c2]);

			// Initial guess of our pressure
			pressure[c] = 0.0f;
		}
	}

	// Jacobi Method Iteration
	for (int k = 0; k < _jNum; k++)
	{
		for (int i = 1; i <= _rows - 2; i++)
		{
			for (int j = 1; j <= _cols - 2; j++)
			{
				// Get grid cell
				int c = idx(i, j);

				// Neighbor Cells
				int c0 = idx(i - 1, j); int c1 = idx(i + 1, j);
				int c2 = idx(i, j - 1); int c3 = idx(i, j + 1);

				// Solve for prressure
				pressure[c] = (pressure[c1] + pressure[c0] + pressure[c3] + pressure[c2] - divergence[c]) / 4.0f;
			}
		}
		setBoundary(pressure, 0);
	}

	//velocity minus grad of Pressure
	for (int i = 1; i <= _rows - 2; i++)
	{
		for (int j = 1; j <= _cols - 2; j++)
		{
			// Get grid cell
			int c = idx(i, j);

			// Neighbor Cells
			int c0 = idx(i - 1, j); int c1 = idx(i + 1, j);
			int c2 = idx(i, j - 1); int c3 = idx(i, j + 1);

			ux[c] -= 0.5f * (pressure[c1] - pressure[c0]);
			uy[c] -= 0.5f * (pressure[c3] - pressure[c2]);
		}
	}
	setBoundary(ux, 1);
	setBoundary(uy, 2);
}

/**
 * resetFields - Resets the velocities/densities
 */
void FluidSolver::resetFields()
{
	// Iterate over the grid
	for (int i = 0; i < _size; i++)
	{
		// Reset Velocity
		ux[i] = 0.0f;
		uy[i] = 0.0f;

		// Reset Density
		d1[i] = 0.0f;
	}
}

/**
 * resetInitialFields - Resets the initial velocities/densities
 */
void FluidSolver::resetInitialFields()
{
	// Iterate over the grid
	for (int i = 0; i < _size; i++)
	{
		// Reset Initial Velocity
		v0x[i] = 0.0f;
		v0y[i] = 0.0f;

		// Reset Initial Density
		d0[i] = 0.0f;
	}
}

void FluidSolver::setBoundary(float *value, int flag)
{
	//for velocity along x-axis
	if (flag == 1)
	{
		for (int i = 1; i <= _rows - 2; i++)
		{
			value[idx(i, 0)] = value[idx(i, 1)];
			value[idx(i, _cols - 1)] = value[idx(i, _cols - 2)];
		}
		for (int j = 1; j <= _cols - 2; j++)
		{
			value[idx(0, j)] = -value[idx(1, j)];
			value[idx(_rows - 1, j)] = -value[idx(_rows - 2, j)];
		}
	}

	//for velocity along y-axis
	if (flag == 2)
	{
		for (int i = 1; i <= _rows - 2; i++)
		{
			value[idx(i, 0)] = -value[idx(i, 1)];
			value[idx(i, _cols - 1)] = -value[idx(i, _cols - 2)];
		}
		for (int j = 1; j <= _cols - 2; j++)
		{
			value[idx(0, j)] = value[idx(1, j)];
			value[idx(_rows - 1, j)] = value[idx(_rows - 2, j)];
		}
	}

	//density
	if (flag == 0)
	{
		for (int i = 1; i <= _rows - 2; i++)
		{
			value[idx(i, 0)] = value[idx(i, 1)];
			value[idx(i, _cols - 1)] = value[idx(i, _cols - 2)];
		}
		for (int j = 1; j <= _cols - 2; j++)
		{
			value[idx(0, j)] = value[idx(1, j)];
			value[idx(_rows - 1, j)] = value[idx(_rows - 2, j)];
		}
	}

	value[idx(0, 0)] = (value[idx(0, 1)] + value[idx(1, 0)]) / 2;
	value[idx(_rows - 1, 0)] = (value[idx(_rows - 2, 0)] + value[idx(_rows - 1, 1)]) / 2;
	value[idx(0, _cols - 1)] = (value[idx(0, _cols - 2)] + value[idx(1, _cols - 1)]) / 2;
	value[idx(_rows - 1, _cols - 1)] = (value[idx(_rows - 2, _cols - 1)] + value[idx(_rows - 1, _cols - 2)]) / 2;
}

/**
 * stepVelocity - moves velocity field forward one timestep
 */
void FluidSolver::stepVelocity()
{
	if (_diffusionRate > 0.0f)
	{
		SWAP(v0x, ux);
		SWAP(v0y, uy);
		diffusion(ux, v0x, _diffusionRate);
		diffusion(uy, v0y, _diffusionRate);
	}

	projection();

	SWAP(v0x, ux);
	SWAP(v0y, uy);
	advection(ux, v0x, v0x, v0y);
	advection(uy, v0y, v0x, v0y);

	projection();
}

/**
 * stepDensity - moves densities forward one timestep
 */
void FluidSolver::stepDensity()
{
	if (_viscosity > 0.0)
	{
		// Diffusion
		SWAP(d0, d1);
		diffusion(d1, d0, _viscosity);
	}

	// Advection
	SWAP(d0, d1);
	advection(d1, d0, ux, uy);
}

/**
 * idx - Gets index into grid
 */
int FluidSolver::idx(int i, int j)
{
	return j * _rows + i;
}

/**
 * setInitVelocity - Sets the Initial Velocity values
 */
void FluidSolver::setInitVelocity(int i, int j, float xVel, float yVel)
{
	// Get grid cell
	int c = idx(i, j);

	// Set Velocity
	v0x[c] = xVel;
	v0y[c] = yVel;
}

/**
 * setInitDensity
 */
void FluidSolver::setInitDensity(int i, int j, float density)
{
	// Get grid cell
	int c = idx(i, j);

	// Set Density
	d0[c] = density;
}

/**
 * setTimestep - Sets the Timestep
 */
void FluidSolver::setTimestep(float dt)
{
	_dt = dt;
}

/**
 * setViscosity - Sets the Viscosity
 */
void FluidSolver::setViscosity(float visc)
{
	_viscosity = visc;
}

/**
 * setDiffusionRate - Sets the Rate of Diffusion for Velocity Field
 */
void FluidSolver::setDiffusionRate(float diffRate)
{
	_diffusionRate = diffRate;
}

/**
 * setJacobiIterations
 */
void FluidSolver::setJacobiIterations(int jNum)
{
	_jNum = jNum;
}

// Getters
point* FluidSolver::getPoints() { return p; }
int FluidSolver::getRows() { return _rows;  }
int FluidSolver::getCols() { return _cols; }
int FluidSolver::getSize() { return _size; }