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

#define SWAP(v0, v1){velocity* tmp = v0; v0 = v1; v1 = tmp;}


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
	void stepVelocity(float visc, float dt);
	void stepDensity();

	void init(int rows, int cols);
	void addSource();
	//void diffusion(velocity* u1, velocity* u0, float visc, float dt);
	void projection();
	void advection(velocity* u1, velocity* u0, float dt);
	void setBoundary();
	int idx(int i, int j);

	// Cleaners
	void reset();
	void resetInitialVelocities();

	// Setters
	void setInitVelocity(int i, int j, float xVel, float yVel);

	// Getters
	point* getPositions();
	velocity* getVelocities();
	int getRows();
	int getCols();
	int getSize();

private:
	// Member Variables
	int rows;
	int cols;
	int size;

	// Grid Points
	point* g;

	// Velocity Arrays
	velocity *u0;
	velocity *u1;

	// Position Array
	point* pos;

	// Time Step
	float t;
};

