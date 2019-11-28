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
	pt = new point[_size];

	// Initialize Velocity Arrays
	vx = new float[_size];
	vy = new float[_size];
	v0x = new float[_size];
	v0y = new float[_size];

	// Initialize Initial Velocity Arrays
	_u0 = new velocity[_size];
	_u1 = new velocity[_size];

	// Initialize Density Arrays
	d1 = new float[_size];
	d0 = new float[_size];

	// Initialize Parameters
	init();
}

/**
 * Destructor
 */
FluidSolver::~FluidSolver()
{
	free(div);
	free(p);
}

/**
 * init
 */
void FluidSolver::init()
{
	h = 1.0f;
	minX = 1.0f;
	maxX = _rows - 1.0f;
	minY = 1.0f;
	maxY = _cols - 1.0f;

	running = 1;

	div = (float *)malloc(sizeof(float)*_size);
	p = (float *)malloc(sizeof(float)*_size);


	for (int i = 0; i < _rows; i++)
	{
		for (int j = 0; j < _cols; j++)
		{
			pt[idx(i, j)].x = (float)i + 0.5f;
			pt[idx(i, j)].y = (float)j + 0.5f;
		}
	}
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
		vx[i] = 0.0f;
		vy[i] = 0.0f;

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
		_u0[i].x = 0.0f;
		_u0[i].y = 0.0f;

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
 * projection
 */
void FluidSolver::projection()
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

			// Calculate Divergence
			div[c] = 0.5f * (vx[c1] - vx[c0] + vy[c3] - vy[c2]);

			// Initial guess of our pressure
			p[c] = 0.0f;
		}
	}
	setBoundary(div, 0);
	setBoundary(p, 0);

	// Jacobi Method Iteration
	for (int k = 0; k < 25; k++)
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
				p[c] = (p[c1] + p[c0] + p[c3] + p[c2] - div[c]) / 4.0f;
			}
		}
		setBoundary(p, 0);
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

			vx[c] -= 0.5f * (p[c1] - p[c0]);
			vy[c] -= 0.5f * (p[c3] - p[c2]);
		}
	}
	setBoundary(vx, 1);
	setBoundary(vy, 2);
}

void FluidSolver::advection(float *value, float *value0, float *u, float *v, int flag)
{
	float oldX;
	float oldY;
	int i0;
	int i1;
	int j0;
	int j1;
	float wL;
	float wR;
	float wB;
	float wT;

	for (int i = 1; i <= _rows - 2; i++)
	{
		for (int j = 1; j <= _cols - 2; j++)
		{
			oldX = pt[idx(i, j)].x - u[idx(i, j)] * _dt;
			oldY = pt[idx(i, j)].y - v[idx(i, j)] * _dt;

			if (oldX < minX) oldX = minX;
			if (oldX > maxX) oldX = maxX;
			if (oldY < minY) oldY = minY;
			if (oldY > maxY) oldY = maxY;

			i0 = (int)(oldX - 0.5f);
			j0 = (int)(oldY - 0.5f);
			i1 = i0 + 1;
			j1 = j0 + 1;

			wL = pt[idx(i1, j0)].x - oldX;
			wR = 1.0f - wL;
			wB = pt[idx(i0, j1)].y - oldY;
			wT = 1.0f - wB;

			value[idx(i, j)] = wB * (wL*value0[idx(i0, j0)] + wR * value0[idx(i1, j0)]) +
				wT * (wL*value0[idx(i0, j1)] + wR * value0[idx(i1, j1)]);
		}
	}

	setBoundary(value, flag);
}

void FluidSolver::diffusion(float *value, float *value0, float rate, int flag)
{
	for (int i = 0; i < _size; i++) value[i] = 0.0f;
	float a = rate * _dt;

	for (int k = 0; k < 20; k++)
	{
		for (int i = 1; i <= _rows - 2; i++)
		{
			for (int j = 1; j <= _cols - 2; j++)
			{
				value[idx(i, j)] = (value0[idx(i, j)] + a * (value[idx(i + 1, j)] + value[idx(i - 1, j)] + value[idx(i, j + 1)] + value[idx(i, j - 1)])) / (4.0f*a + 1.0f);
			}
		}
		setBoundary(value, flag);
	}
}

/**
 * diffuseVelocity
 */
void FluidSolver::diffuseVelocity(velocity *u1, velocity *u0, float rate, int flag)
{
	for (int i = 0; i < _size; i++)
	{
		u1[i].x = 0.0f;
		u1[i].y = 0.0f;
	}

	float a = rate * _dt;

	for (int k = 0; k < 20; k++)
	{
		for (int i = 1; i <= _rows - 2; i++)
		{
			for (int j = 1; j <= _cols - 2; j++)
			{
				// Get grid cell
				int c = idx(i, j);

				// Neighbor Cells
				int c0 = idx(i - 1, j);
				int c1 = idx(i + 1, j);
				int c2 = idx(i, j - 1);
				int c3 = idx(i, j + 1);

				u1[c].x = (u0[c].x + a * (u1[c1].x + u1[c0].x + u1[c3].x + u1[c2].x)) / (4.0f*a + 1.0f);
				u1[c].y = (u0[c].y + a * (u1[c1].y + u1[c0].y + u1[c3].y + u1[c2].y)) / (4.0f*a + 1.0f);
			}
		}
		//setBoundary(value, flag);
	}
}

/**
 * addForce - Adds the Initial Force to velocity/density
 */
void FluidSolver::addForce()
{
	// Iterate over the Grid
	for (int i = 1; i <= _size; i++)
	{
		// Set Velocity
		vx[i] += v0x[i];
		vy[i] += v0y[i];
		_u1[i].x += _u0[i].x;
		_u1[i].y += _u0[i].y;

		// Set Density
		d1[i] += d0[i];
	}
}

/**
 * stepVelocity - moves velocity field forward one timestep
 */
void FluidSolver::stepVelocity()
{
	if (_diffusionRate > 0.0f)
	{
		SWAP(v0x, vx);
		SWAP(v0y, vy);
		diffusion(vx, v0x, _diffusionRate, 1);
		diffusion(vy, v0y, _diffusionRate, 2);
	}

	projection();

	SWAP(v0x, vx);
	SWAP(v0y, vy);
	advection(vx, v0x, v0x, v0y, 1);
	advection(vy, v0y, v0x, v0y, 2);

	projection();
}

/**
 * stepDensity - moves densities forward one timestep
 */
void FluidSolver::stepDensity()
{
	// Advection
	SWAP(d0, d1);
	advection(d1, d0, vx, vy, 0);

	// Diffusion
	SWAP(d0, d1);
	diffusion(d1, d0, _viscosity, 0);
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
	_u0[c].x = xVel;
	_u0[c].y = yVel;
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

// Getters
point* FluidSolver::getPoints() { return pt; }
int FluidSolver::getRows() { return _rows;  }
int FluidSolver::getCols() { return _cols; }
int FluidSolver::getSize() { return _size; }