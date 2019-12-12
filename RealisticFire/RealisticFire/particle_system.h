#pragma once

#include <vector>
#include "particle.h"

class particle_system
{
	vector<particle> particles;
	vec3d gravity_point;

public:
	//construct system given n number of particles
	particle_system(int);

	//Function to advance state of particle system by time t in ms
	void advance(float);

	//Function to set gravity point
	void set_gravity(vec3d = vec3d(0, 50, 0));

	//Function to add particles
	bool add_particles(int);

	//Function to delete particles
	bool delete_particles(int);

	bool delete_particle();

	bool add_particle(int);

	//Function to draw particles
	void draw(int tex);

	~particle_system(void);
};

