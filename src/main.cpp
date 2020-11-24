#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <random>


// This Struct represents a point particle in 3D space
struct particle {
	double x;
	double y;
	double z;

	particle(double x, double y, double z)
		: x(x), y(y), z(z) {};
};

// Distance between two particles
double distance(particle p1, particle p2) {
	return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));
}


/*
	Objects of this class represent instances of an active molecular dynamics simulation

		
	@param n_particles (int): Number of particles in the simulation
	@param box_dim (double): Dimensions of the simulation's 3D box
	@param timestep (double): Timestep of the simulation
	@param initial_velocity (double): Max initial velocity of the particles per timestep
	@param sigma (double): Size of the particles in distance units
	@param mass (double): Mass of the particles
*/
class simulation {
private:
	// Simulation parameters
	int n_particles;
	double box_x, box_y, box_z;
	double timestep;
	int n_steps;
	
	// Particle parameters
	double initial_velocity;
	double sigma;	// Size of particles (distance units)
	double mass;

	// Simulation Particles
	particle* particles_t0;
	particle* particles_t1;
	particle* particles_t2;

	void particle_init() {

		bool colision;

		// Random Number from Hardware
		std::random_device dev;
		// Seed Pesudo-Random Generator
		std::mt19937 gen(dev());
		// Real distribution from 0 to 1
		std::uniform_real_distribution<> distr(0,1);

		// Initialize 1st Particle
		particles_t0[0] = particle(
			distr(gen) * box_x, 
			distr(gen) * box_y, 
			distr(gen) * box_z);

		particles_t1[0] = particle(
			initial_velocity * (1. - 2. * distr(gen)) + particles_t0[0].x, 
			initial_velocity * (1. - 2. * distr(gen)) + particles_t0[0].y, 
			initial_velocity * (1. - 2. * distr(gen)) + particles_t0[0].z);

		for (int i = 1; i < this->n_particles; i++) {

			colision = true;
			 
			while (colision) {

				colision = false;

				// Generate Particle
				particles_t0[i] = particle(
					distr(gen) * box_x,
					distr(gen) * box_y,
					distr(gen) * box_z);

				// Check if particle is to close to another
				for (int j = 0; j < i; j++) {
					if (this->nearest_image(particles_t0[i], particles_t0[j]) < 1.5 * this->sigma) {
						colision = true;
					}
				}
			}

			// Give random velocity to the particle
			particles_t1[i] = particle(
				initial_velocity * (1. - 2. * distr(gen)) + particles_t0[i].x,
				initial_velocity * (1. - 2. * distr(gen)) + particles_t0[i].y,
				initial_velocity * (1. - 2. * distr(gen)) + particles_t0[i].z);

		}
	}

public:

	// Initialize Simulation Instance
	simulation(int n_particles, double box_dim, double timestep, int n_step, double initial_velocity, double sigma, double mass)
		: n_particles(n_particles), box_x(box_dim), box_y(box_dim), box_z(box_dim), n_steps(n_step), initial_velocity(initial_velocity), sigma(sigma), mass(mass)
	{
		// Particle Array Initialization for 3 Time Instances
		
		// Time Instance 0
		particles_t0 = (particle*) malloc(this->n_particles * sizeof(particle));
		// Time Instance 1
		particles_t1 = (particle*) malloc(this->n_particles * sizeof(particle));
		// Time Instance 2
		particles_t2 = (particle*) malloc(this->n_particles * sizeof(particle));

		// Create Particles
		this->particle_init();
	}

	/**
	 * @brief Calculates the nearest distance between two particles in a 3D periodic system.
	 *
	 * The function implements the nearest image convention (NIC), commonly used in the simulation of
	 * periodic systems such as molecular dynamics simulations. It computes the shortest distance
	 * between two points in a system with periodic boundary conditions (a point exiting one side of the
	 * box re-enters from the opposite side).
	 *
	 * @param p1 The first particle, containing x, y, and z coordinates.
	 * @param p2 The second particle, also containing x, y, and z coordinates.
	 * @return The shortest distance between p1 and p2 considering periodic boundary conditions.
	 */	
	double nearest_image(particle p1, particle p2) {

		double x, y, z;

		x = p1.x - p2.x;
		x = x - this->box_x * round(x / this->box_x);
		y = p1.y - p2.y;
		y = y - this->box_y * round(y / this->box_y);
		z = p1.z - p2.z;
		z = z - this->box_z * round(z / this->box_z);

		return sqrt(x * x + y * y + z * z);
	}

};


int main(int argc, int** argv) {
	return 0;
}