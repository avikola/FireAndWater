// Standard Libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

// Headers
#include "Fluid.h"

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
	_u0 = new velocity[_size];
	_u1 = new velocity[_size];

	// Initialize Density Arrays
	_d1 = new float[_size];
	_d0 = new float[_size];

	// Initialize Divergence
	divergence = new float[_size];

	// Initialize Pressure
	_pressure = new float[_size];

	minX = 1.0f;
	maxX = _rows - 1.0f;
	minY = 1.0f;
	maxY = _cols - 1.0f;

	// Iterate over rows
	for (int i = 0; i < _rows; ++i)
	{
		// Iterate over cols
		for (int j = 0; j < _cols; ++j)
		{
			// Get grid cell
			int c = idx(i, j);

			p[c].x = (float)i + 0.5f;
			p[c].y = (float)j + 0.5f;
		}
	}
}

/**
 * Destructor
 */
FluidSolver::~FluidSolver() {}

/**
 * addForce - Adds the Initial Force to velocity/density
 */
void FluidSolver::addForce()
{
	// Iterate over the Grid
	for (int i = 1; i <= _size; ++i)
	{
		// Set Velocity
		_u1[i].x += _u0[i].x;
		_u1[i].y += _u0[i].y;

		// Set Density
		_d1[i] += _d0[i];
	}
}


/**
 * advectDensity
 */
void FluidSolver::advectDensity(float* d1, float* d0, velocity* u1)
{
	for (int i = 0; i < _rows; ++i)
	{
		for (int j = 0; j < _cols; ++j)
		{
			// Get grid cell
			int c = idx(i, j);

			float oldX = p[c].x - u1[c].x * _dt;
			float oldY = p[c].y - u1[c].y * _dt;

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

			d1[idx(i, j)] = wB * (wL*d0[idx(i0, j0)] + wR * d0[idx(i1, j0)]) +
							wT * (wL*d0[idx(i0, j1)] + wR * d0[idx(i1, j1)]);
		}
	}
}

/**
 * advectVelocity
 */
void FluidSolver::advectVelocity(velocity* u1, velocity* u0)
{
	for (int i = 0; i < _rows; ++i)
	{
		for (int j = 0; j < _cols; ++j)
		{
			// Get grid cell
			int c = idx(i, j);

			float oldX = p[c].x - u0[c].x * _dt;
			float oldY = p[c].y - u0[c].y * _dt;

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

			u1[idx(i, j)].x = wB * (wL*u0[idx(i0, j0)].x + wR * u0[idx(i1, j0)].x) +
				wT * (wL*u0[idx(i0, j1)].x + wR * u0[idx(i1, j1)].x);

			u1[idx(i, j)].y = wB * (wL*u0[idx(i0, j0)].y + wR * u0[idx(i1, j0)].y) +
				wT * (wL*u0[idx(i0, j1)].y + wR * u0[idx(i1, j1)].y);
		}
	}
}

/**
 * diffusion
 */
void FluidSolver::diffusion(float *v, float *v0, float rate)
{
	for (int i = 0; i < _size; ++i)
		v[i] = 0.0f;
	
	float a = rate * _dt;

	// Jacobi Method Iteration
	for (int k = 0; k < _jNum; ++k)
	{
		// Iterate over rows in grid
		for (int i = 0; i < _rows; ++i)
		{
			// Iterate over cols in grid
			for (int j = 0; j < _cols; ++j)
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
 * diffuseVelocity
 */
void FluidSolver::diffuseVelocity(velocity* u1, velocity* u0, float rate)
{
	for (int i = 0; i < _size; ++i)
	{
		u1[i].x = 0.0f;
		u1[i].y = 0.0f;
	}
	float a = rate * _dt;

	// Jacobi Method Iteration
	for (int k = 0; k < _jNum; ++k)
	{
		// Iterate over rows in grid
		for (int i = 0; i < _rows; ++i)
		{
			// Iterate over cols in grid
			for (int j = 0; j < _cols; ++j)
			{
				// Get grid cell
				int c = idx(i, j);

				// Check for Boundary cases
				int i0 = (_rows + ((i - 1) % _rows)) + _rows;
				int j0 = (_cols + ((j - 1) % _cols)) + _cols;
				int i1 = (i + 1) % _rows;
				int j1 = (j + 1) % _cols;

				// Neighbor grid Cells
				int c0 = idx(i0, j); int c1 = idx(i1, j);
				int c2 = idx(i, j0); int c3 = idx(i, j1);

				u1[c].x = (u0[c].x + a * (u1[c1].x + u1[c0].x + u1[c3].x + u1[c2].x)) / (4.0f*a + 1.0f);
				u1[c].y = (u0[c].y + a * (u1[c1].y + u1[c0].y + u1[c3].y + u1[c2].y)) / (4.0f*a + 1.0f);
			}
		}
	}
}

/**
 * projection
 */
void FluidSolver::projection()
{
	// Iterate over rows
	for (int i = 0; i < _rows; ++i)
	{
		// Iterate over cols
		for (int j = 0; j < _cols; ++j)
		{
			// Get grid cell
			int c = idx(i, j);

			// Neighbor Cells
			int c0 = idx(i - 1, j); int c1 = idx(i + 1, j);
			int c2 = idx(i, j - 1); int c3 = idx(i, j + 1);

			// Calculate Divergence
			divergence[c] = 0.5f * (_u1[c1].x - _u1[c0].x + _u1[c3].y - _u1[c2].y);

			// Initial Pressure
			_pressure[c] = 0.0f;
		}
	}

	// Calculate Pressure
	computePressure();

	// Iterate over rows
	for (int i = 1; i <= _rows - 2; ++i)
	{
		// Iterate over cols
		for (int j = 1; j <= _cols - 2; ++j)
		{
			// Get grid cell
			int c = idx(i, j);

			// Neighbor Cells
			int c0 = idx(i - 1, j); int c1 = idx(i + 1, j);
			int c2 = idx(i, j - 1); int c3 = idx(i, j + 1);

			// Velocity minus gradiant of Pressure
			_u1[c].x -= 0.5f * (_pressure[c1] - _pressure[c0]);
			_u1[c].y -= 0.5f * (_pressure[c3] - _pressure[c2]);
		}
	}

	setVelocityBoundary(_u1);
}

/**
 * computePressure - Calculates the Pressure
 */
void FluidSolver::computePressure()
{
	// Jacobi Method Iteration
	for (int k = 0; k < _jNum; ++k)
	{
		// Iterate over rows
		for (int i = 1; i <= _rows - 2; ++i)
		{
			// Iterate over cols
			for (int j = 1; j <= _cols - 2; ++j)
			{
				// Get grid cell
				int c = idx(i, j);

				// Neighbor Cells
				int c0 = idx(i - 1, j); int c1 = idx(i + 1, j);
				int c2 = idx(i, j - 1); int c3 = idx(i, j + 1);

				// Iteratively solve for pressure
				_pressure[c] = (_pressure[c1] + _pressure[c0] + _pressure[c3] + _pressure[c2] - divergence[c]) / 4.0f;
			}
		}

		// Set Pressure Boundary
		setPressureBoundary(_pressure);
	}
}

/**
 * resetFields - Resets the velocities/densities
 */
void FluidSolver::resetFields()
{
	// Iterate over the grid
	for (int i = 0; i < _size; ++i)
	{
		// Reset Velocity
		_u1[i].x = 0.0f;
		_u1[i].y = 0.0f;

		// Reset Density
		_d1[i] = 0.0f;
	}
}

/**
 * resetInitialFields - Resets the initial velocities/densities
 */
void FluidSolver::resetInitialFields()
{
	// Iterate over the grid
	for (int i = 0; i < _size; ++i)
	{
		// Reset Initial Velocity
		_u0[i].x = 0.0f;
		_u0[i].y = 0.0f;

		// Reset Initial Density
		_d0[i] = 0.0f;
	}
}

/**
 * setPressureBoundary
 */
void FluidSolver::setPressureBoundary(float *pressure)
{
	for (int i = 1; i <= _rows - 2; ++i)
	{
		pressure[idx(i, 0)] = pressure[idx(i, 1)];
		pressure[idx(i, _cols - 1)] = pressure[idx(i, _cols - 2)];
	}

	for (int j = 1; j <= _cols - 2; ++j)
	{
		pressure[idx(0, j)] = pressure[idx(1, j)];
		pressure[idx(_rows - 1, j)] = pressure[idx(_rows - 2, j)];
	}

	// Update Pressure at Corners
	pressure[idx(0, 0)] = (pressure[idx(0, 1)] + pressure[idx(1, 0)]) / 2;
	pressure[idx(_rows - 1, 0)] = (pressure[idx(_rows - 2, 0)] + pressure[idx(_rows - 1, 1)]) / 2;
	pressure[idx(0, _cols - 1)] = (pressure[idx(0, _cols - 2)] + pressure[idx(1, _cols - 1)]) / 2;
	pressure[idx(_rows - 1, _cols - 1)] = (pressure[idx(_rows - 2, _cols - 1)] + pressure[idx(_rows - 1, _cols - 2)]) / 2;
}

/**
 * setVelocityBoundary
 */
void FluidSolver::setVelocityBoundary(velocity* u1)
{
	for (int i = 1; i < _rows - 1; ++i)
	{
		// Update X-Axis Velocity
		u1[idx(i, 0)].x = u1[idx(i, 1)].x;
		u1[idx(i, _cols - 1)].x = u1[idx(i, _cols - 2)].x;

		// Update Y-Axis Velocity
		u1[idx(i, 0)].y = -u1[idx(i, 1)].y;
		u1[idx(i, _cols - 1)].y = -u1[idx(i, _cols - 2)].y;
	}

	for (int j = 1; j < _cols - 1; ++j)
	{
		// Update X-Axis Velocity
		u1[idx(0, j)].x = -u1[idx(1, j)].x;
		u1[idx(_rows - 1, j)].x = -u1[idx(_rows - 2, j)].x;

		// Update Y-Axis Velocity
		u1[idx(0, j)].y = u1[idx(1, j)].y;
		u1[idx(_rows - 1, j)].y = u1[idx(_rows - 2, j)].y;
	}

	// Update Corners X Velocity
	u1[idx(0, 0)].x = (u1[idx(0, 1)].x + u1[idx(1, 0)].x) / 2;
	u1[idx(_rows - 1, 0)].x = (u1[idx(_rows - 2, 0)].x + u1[idx(_rows - 1, 1)].x) / 2;
	u1[idx(0, _cols - 1)].x = (u1[idx(0, _cols - 2)].x + u1[idx(1, _cols - 1)].x) / 2;
	u1[idx(_rows - 1, _cols - 1)].x = (u1[idx(_rows - 2, _cols - 1)].x + u1[idx(_rows - 1, _cols - 2)].x) / 2;

	// Update Corners Y Velocity
	u1[idx(0, 0)].y = (u1[idx(0, 1)].y + u1[idx(1, 0)].y) / 2;
	u1[idx(_rows - 1, 0)].y = (u1[idx(_rows - 2, 0)].y + u1[idx(_rows - 1, 1)].y) / 2;
	u1[idx(0, _cols - 1)].y = (u1[idx(0, _cols - 2)].y + u1[idx(1, _cols - 1)].y) / 2;
	u1[idx(_rows - 1, _cols - 1)].y = (u1[idx(_rows - 2, _cols - 1)].y + u1[idx(_rows - 1, _cols - 2)].y) / 2;
}

/**
 * stepVelocity - moves velocity field forward one timestep
 */
void FluidSolver::stepVelocity()
{
	// Diffusion
	if (_diffusionRate > 0.0f)
	{
		swapVelocity(_u1, _u0);
		diffuseVelocity(_u1, _u0, _diffusionRate);
	}

	// Projection
	projection();

	// Advection
	swapVelocity(_u1, _u0);
	advectVelocity(_u1, _u0);
}

/**
 * stepDensity - moves densities forward one timestep
 */
void FluidSolver::stepDensity()
{
	if (_viscosity > 0.0)
	{
		// Diffusion
		SWAP(_d0, _d1);
		diffusion(_d1, _d0, _viscosity);
	}

	// Advection
	SWAP(_d0, _d1);
	advectDensity(_d1, _d0, _u1);
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
	_u0[c].x = xVel;
	_u0[c].y = yVel;
}

/**
 * swapVelocity
 */
void FluidSolver::swapVelocity(velocity* u1, velocity* u0)
{
	for (int i = 1; i < _size; ++i)
	{
		float tmpX = u0[i].x;
		float tmpY = u0[i].y;

		u0[i].x = u1[i].x;
		u0[i].y = u1[i].y;

		u1[i].x = tmpX;
		u1[i].y = tmpY;
	}
}

/**
 * setInitDensity
 */
void FluidSolver::setInitDensity(int i, int j, float density)
{
	// Get grid cell
	int c = idx(i, j);

	// Set Density
	_d0[c] = density;
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
velocity* FluidSolver::getVelocities() { return _u1; }
float* FluidSolver::getDensities() { return _d1; }
int FluidSolver::getRows() { return _rows;  }
int FluidSolver::getCols() { return _cols; }
int FluidSolver::getSize() { return _size; }

float FluidSolver::getDensity(int i, int j, float factor)	
{
	 return (_d1[idx(i - 1, j - 1)] + _d1[idx(i, j - 1)] + _d1[idx(i - 1, j)] + _d1[idx(i, j)]) / factor;
}