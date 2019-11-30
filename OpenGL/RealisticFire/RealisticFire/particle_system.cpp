#include "particle_system.h"
#include <GL/glut.h>
#include <math.h>

particle_system::particle_system(int n)
{
	if (n > MAX_PARTICLES)
		n = MAX_PARTICLES;

	for (int i = 0; i < n; i++)
	{
		particle temp;
		particles.push_back(temp);
	}
}


//Function to advance state of particle system by time t in ms
void particle_system::advance(float time)
{
	vector<particle>::iterator it;
	for (it = particles.begin(); it != particles.end(); it++)
	{
		vec3d force = vec3d((gravity_point - it->get_pos()).unit(), FORCE_MAG);
		it->advance(time, force);
	}
}

//Function to set gravity point
void particle_system::set_gravity(vec3d gravity) {
	gravity_point = gravity;
}

//Function to add particles
bool particle_system::add_particles(int num)
{
	int i;
	for (i = 0; i < num && particles.size() < MAX_PARTICLES; i++)
	{
		particle p;
		particles.push_back(p);
	}
	return (i >= num);
}

//Function to delete particles(least is 0)
bool particle_system::delete_particles(int num)
{
	int i;
	for (i = 0; i < num && particles.size() > 0; i++)
	{
		particles.pop_back();
	}

	return (i >= num);
}

//Function to delete particles(least is 0)
bool particle_system::delete_particle()
{
	int size = particles.size();
	for (int i = 0; i < size; i++)
	{
		particle p;
		vector<particle>::iterator it = particles.begin() + i;
		if (it->get_life() == 0) {
			particles.erase(it);
			particles.push_back(p);
		}
	}

	return true;
}

//Function to draw a particle
void particle_system::draw()
{
	vector<particle>::iterator it;
	for (it = particles.begin(); it != particles.end(); it++) {
		vec3d pos = it->get_pos();
		float k = abs(pos.x) / (.5*LENGTH);
		float k2 = abs(pos.y) / (1.5*LENGTH);
		glColor4f(1, 1 - k, .3 - .3*k2, 1);
		glBegin(GL_TRIANGLES);
		glVertex3f(pos.x - .5 + LENGTH, pos.y - .5 + LENGTH, pos.z);
		glVertex3f(pos.x + .5 + LENGTH, pos.y - .5 + LENGTH, pos.z);
		glVertex3f(pos.x + LENGTH, pos.y + 1 + LENGTH, pos.z);
		glEnd();
	}
}


particle_system::~particle_system(void)
{
}