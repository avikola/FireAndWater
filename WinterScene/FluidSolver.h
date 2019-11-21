
struct point
{
	float x;
	float y;
};

struct velocity
{
	float x;
	float y;
};


class FluidSolver
{
public:
	FluidSolver();
	~FluidSolver();

	enum BoundaryType
	{
		X_VELOCITY,
		Y_VELOCITY,
		DENSITY
	};



	// Step Functions for Velocity/Diffusion
	void stepVelocity();
	void stepDensity();

	void init(int rows, int cols);
	void diffusion(float *d, float *d0, float dt, float h, float r);
	void projection();
	void advection(float *val0, float *val, float *u, float *v, float dt);
	void setBoundary();
	int cellIndex(int i, int j);

	// Getters
	float* getVelocityX();
	float* getVelocityY();
	point* getPositions();

private:
	// Member Variables
	int rows;
	int cols;
	int size;

	// Characteristic
	float *px;
	float *py;

	// Previous X/Y Velocities
	float *v0X;
	float *v0Y;

	// X/Y Velocities
	float *vX;
	float *vY;

	// Position Array
	point *pos;

	// Time Step
	float t;
};

