#include <gl/glut.h>

#include "particle.h"

#define PI 3.14159

#include "math.h"


particle::particle(float _mass, vec3d _velocity, vec3d _position) : mass(_mass), velocity(_velocity), position(_position)
{
	float random = rand_float();
	float random2 = rand_float();

	position.x = cos(2 * PI * random) * 100 * random2;
	position.y = sin(2 * PI * random) * 5 * random2 - 80;

	life = (200 - .1*pow(abs(position.x),2)) + rand_float() * 50;
}

// Function to advance state by time t in ms
void particle::advance(float t, vec3d force)
{
	// Calculating acceleration
	life = life - t;

	if (life < 0) 
		life = 0;
	
	vec3d acc = force / mass;

	// Calculating velocity
	velocity = velocity + acc * (t / 1000.0);

	if (velocity.mag() >= MAX_VELOCITY)
		velocity = vec3d(velocity.unit(), MAX_VELOCITY);

	// Changing position
	position = position + velocity * (t / 500.0);

	if (position.x <= -LENGTH)
		position.x = LENGTH;
	else if (position.x >= LENGTH)
		position.x = -LENGTH;

	if (position.y <= -LENGTH)
		position.y = LENGTH;
	else if (position.y >= LENGTH)
		position.y = -LENGTH;

	if (position.z <= -LENGTH)
		position.z = LENGTH;
	else if (position.z >= LENGTH)
		position.z = -LENGTH;
}

particle::~particle(void) {}

// Function to get position
vec3d particle::get_pos() { return position; }

float particle::get_life() { return life; }
