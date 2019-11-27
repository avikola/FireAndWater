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
 * rows - Rows in Grid
 * cols - Cols in Grid
 * dt   - TimeStep
 * visc - Viscosity
 */
FluidSolver::FluidSolver(int rows, int cols, float dt, float visc)
{
	// Set Dimension Sizes
	_rows = rows;
	_cols = cols;

	// Set Grid Size
	_size = _rows * _cols;

	// Initialize TimeStep
	_dt = dt;

	// Initialize Viscosity
	_viscosity = visc;

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

	//vorticity confinement
	free(vort);
	free(absVort);
	free(gradVortX);
	free(gradVortY);
	free(lenGrad);
	free(vcfx);
	free(vcfy);
}

/**
 * init
 */
void FluidSolver::init()
{
	h = 1.0f;
	simSizeX = (float)_rows;
	simSizeY = (float)_cols;
	minX = 1.0f;
	maxX = _rows - 1.0f;
	minY = 1.0f;
	maxY = _cols - 1.0f;

	running = 1;

	diff = 0.0f;
	vorticity = 0.0f;

	div = (float *)malloc(sizeof(float)*_size);
	p = (float *)malloc(sizeof(float)*_size);

	//vorticity confinement
	vort = (float *)malloc(sizeof(float)*_size);
	absVort = (float *)malloc(sizeof(float)*_size);
	gradVortX = (float *)malloc(sizeof(float)*_size);
	gradVortY = (float *)malloc(sizeof(float)*_size);
	lenGrad = (float *)malloc(sizeof(float)*_size);
	vcfx = (float *)malloc(sizeof(float)*_size);
	vcfy = (float *)malloc(sizeof(float)*_size);

	for (int i = 0; i < _rows; i++)
	{
		for (int j = 0; j < _cols; j++)
		{
			pt[cIdx(i, j)].x = (float)i + 0.5f;
			pt[cIdx(i, j)].y = (float)j + 0.5f;
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
			value[cIdx(i, 0)] = value[cIdx(i, 1)];
			value[cIdx(i, _cols - 1)] = value[cIdx(i, _cols - 2)];
		}
		for (int j = 1; j <= _cols - 2; j++)
		{
			value[cIdx(0, j)] = -value[cIdx(1, j)];
			value[cIdx(_rows - 1, j)] = -value[cIdx(_rows - 2, j)];
		}
	}

	//for velocity along y-axis
	if (flag == 2)
	{
		for (int i = 1; i <= _rows - 2; i++)
		{
			value[cIdx(i, 0)] = -value[cIdx(i, 1)];
			value[cIdx(i, _cols - 1)] = -value[cIdx(i, _cols - 2)];
		}
		for (int j = 1; j <= _cols - 2; j++)
		{
			value[cIdx(0, j)] = value[cIdx(1, j)];
			value[cIdx(_rows - 1, j)] = value[cIdx(_rows - 2, j)];
		}
	}

	//density
	if (flag == 0)
	{
		for (int i = 1; i <= _rows - 2; i++)
		{
			value[cIdx(i, 0)] = value[cIdx(i, 1)];
			value[cIdx(i, _cols - 1)] = value[cIdx(i, _cols - 2)];
		}
		for (int j = 1; j <= _cols - 2; j++)
		{
			value[cIdx(0, j)] = value[cIdx(1, j)];
			value[cIdx(_rows - 1, j)] = value[cIdx(_rows - 2, j)];
		}
	}

	value[cIdx(0, 0)] = (value[cIdx(0, 1)] + value[cIdx(1, 0)]) / 2;
	value[cIdx(_rows - 1, 0)] = (value[cIdx(_rows - 2, 0)] + value[cIdx(_rows - 1, 1)]) / 2;
	value[cIdx(0, _cols - 1)] = (value[cIdx(0, _cols - 2)] + value[cIdx(1, _cols - 1)]) / 2;
	value[cIdx(_rows - 1, _cols - 1)] = (value[cIdx(_rows - 2, _cols - 1)] + value[cIdx(_rows - 1, _cols - 2)]) / 2;
}

void FluidSolver::projection()
{
	for (int i = 1; i <= _rows - 2; i++)
	{
		for (int j = 1; j <= _cols - 2; j++)
		{
			div[cIdx(i, j)] = 0.5f * (vx[cIdx(i + 1, j)] - vx[cIdx(i - 1, j)] + vy[cIdx(i, j + 1)] - vy[cIdx(i, j - 1)]);
			p[cIdx(i, j)] = 0.0f;;
		}
	}
	setBoundary(div, 0);
	setBoundary(p, 0);

	//projection iteration
	for (int k = 0; k < 20; k++)
	{
		for (int i = 1; i <= _rows - 2; i++)
		{
			for (int j = 1; j <= _cols - 2; j++)
			{
				p[cIdx(i, j)] = (p[cIdx(i + 1, j)] + p[cIdx(i - 1, j)] + p[cIdx(i, j + 1)] + p[cIdx(i, j - 1)] - div[cIdx(i, j)]) / 4.0f;
			}
		}
		setBoundary(p, 0);
	}

	//velocity minus grad of Pressure
	for (int i = 1; i <= _rows - 2; i++)
	{
		for (int j = 1; j <= _cols - 2; j++)
		{
			vx[cIdx(i, j)] -= 0.5f*(p[cIdx(i + 1, j)] - p[cIdx(i - 1, j)]);
			vy[cIdx(i, j)] -= 0.5f*(p[cIdx(i, j + 1)] - p[cIdx(i, j - 1)]);
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
			oldX = pt[cIdx(i, j)].x - u[cIdx(i, j)] * _dt;
			oldY = pt[cIdx(i, j)].y - v[cIdx(i, j)] * _dt;

			if (oldX < minX) oldX = minX;
			if (oldX > maxX) oldX = maxX;
			if (oldY < minY) oldY = minY;
			if (oldY > maxY) oldY = maxY;

			i0 = (int)(oldX - 0.5f);
			j0 = (int)(oldY - 0.5f);
			i1 = i0 + 1;
			j1 = j0 + 1;

			wL = pt[cIdx(i1, j0)].x - oldX;
			wR = 1.0f - wL;
			wB = pt[cIdx(i0, j1)].y - oldY;
			wT = 1.0f - wB;

			value[cIdx(i, j)] = wB * (wL*value0[cIdx(i0, j0)] + wR * value0[cIdx(i1, j0)]) +
				wT * (wL*value0[cIdx(i0, j1)] + wR * value0[cIdx(i1, j1)]);
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
				value[cIdx(i, j)] = (value0[cIdx(i, j)] + a * (value[cIdx(i + 1, j)] + value[cIdx(i - 1, j)] + value[cIdx(i, j + 1)] + value[cIdx(i, j - 1)])) / (4.0f*a + 1.0f);
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
				int c = cIdx(i, j);

				// Neighbor Cells
				int c0 = cIdx(i - 1, j);
				int c1 = cIdx(i + 1, j);
				int c2 = cIdx(i, j - 1);
				int c3 = cIdx(i, j + 1);

				u1[c].x = (u0[c].x + a * (u1[c1].x + u1[c0].x + u1[c3].x + u1[c2].x)) / (4.0f*a + 1.0f);
				u1[c].y = (u0[c].y + a * (u1[c1].y + u1[c0].y + u1[c3].y + u1[c2].y)) / (4.0f*a + 1.0f);
			}
		}
		//setBoundary(value, flag);
	}
}

void FluidSolver::vortConfinement()
{
	for (int i = 1; i <= _rows - 2; i++)
	{
		for (int j = 1; j <= _cols - 2; j++)
		{
			vort[cIdx(i, j)] = 0.5f*(vy[cIdx(i + 1, j)] - vy[cIdx(i - 1, j)] - vx[cIdx(i, j + 1)] + vx[cIdx(i, j - 1)]);
			if (vort[cIdx(i, j)] >= 0.0f) absVort[cIdx(i, j)] = vort[cIdx(i, j)];
			else absVort[cIdx(i, j)] = -vort[cIdx(i, j)];
		}
	}
	setBoundary(vort, 0);
	setBoundary(absVort, 0);

	for (int i = 1; i <= _rows - 2; i++)
	{
		for (int j = 1; j <= _cols - 2; j++)
		{
			gradVortX[cIdx(i, j)] = 0.5f*(absVort[cIdx(i + 1, j)] - absVort[cIdx(i - 1, j)]);
			gradVortY[cIdx(i, j)] = 0.5f*(absVort[cIdx(i, j + 1)] - absVort[cIdx(i, j - 1)]);
			lenGrad[cIdx(i, j)] = sqrt(gradVortX[cIdx(i, j)] * gradVortX[cIdx(i, j)] + gradVortY[cIdx(i, j)] * gradVortY[cIdx(i, j)]);
			if (lenGrad[cIdx(i, j)] < 0.01f)
			{
				vcfx[cIdx(i, j)] = 0.0f;
				vcfy[cIdx(i, j)] = 0.0f;
			}
			else
			{
				vcfx[cIdx(i, j)] = gradVortX[cIdx(i, j)] / lenGrad[cIdx(i, j)];
				vcfy[cIdx(i, j)] = gradVortY[cIdx(i, j)] / lenGrad[cIdx(i, j)];
			}
		}
	}
	setBoundary(vcfx, 0);
	setBoundary(vcfy, 0);

	for (int i = 1; i <= _rows - 2; i++)
	{
		for (int j = 1; j <= _cols - 2; j++)
		{
			vx[cIdx(i, j)] += vorticity * (vcfy[cIdx(i, j)] * vort[cIdx(i, j)]);
			vy[cIdx(i, j)] += vorticity * (-vcfx[cIdx(i, j)] * vort[cIdx(i, j)]);
		}
	}

	setBoundary(vx, 1);
	setBoundary(vy, 2);
}

/**
 * addSource
 */
void FluidSolver::addSource()
{
	for (int i = 1; i <= _rows - 2; i++)
	{
		for (int j = 1; j <= _cols - 2; j++)
		{
			// Get grid cell
			int c = cIdx(i, j);

			// Set Velocity
			vx[c] += v0x[c];
			vy[c] += v0y[c];
			_u1[c].x += _u0[c].x;
			_u1[c].y += _u0[c].y;


			// Set Density
			d1[c] += d0[c];
		}
	}

	setBoundary(vx, 1);
	setBoundary(vy, 2);
	setBoundary(d1, 0);
}

/**
 * stepVelocity - moves velocity field forward one timestep
 */
void FluidSolver::stepVelocity()
{
	if (diff > 0.0f)
	{
		SWAP(v0x, vx);
		SWAP(v0y, vy);
		diffusion(vx, v0x, diff, 1);
		diffusion(vy, v0y, diff, 2);
	}

	projection();

	SWAP(v0x, vx);
	SWAP(v0y, vy);
	advection(vx, v0x, v0x, v0y, 1);
	advection(vy, v0y, v0x, v0y, 2);

	projection();

	vortConfinement();
}

/**
 * stepDensity - moves densities forward one timestep
 */
void FluidSolver::stepDensity()
{
	if (_viscosity > 0.0f)
	{
		SWAP(d0, d1);
		diffusion(d1, d0, _viscosity, 0);
	}
	SWAP(d0, d1);
	advection(d1, d0, vx, vy, 0);
}

/**
 * setInitVelocity - Sets the Initial Velocity values
 */
void FluidSolver::setInitVelocity(int i, int j, float xVel, float yVel)
{
	// Get grid cell
	int c = cIdx(i, j);

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
	int c = cIdx(i, j);

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

// Getters
point* FluidSolver::getPoints() { return pt; }
int FluidSolver::getRows() { return _rows;  }