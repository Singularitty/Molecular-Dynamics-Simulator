#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include "particle.hpp"
#include "../Potentials/PotentialField.hpp"
#include "../Output/SimulationOutput.hpp"

/*
	Objects of this class represent instances of a molecular dynamics simulation

	@param n_particles (int): Number of particles in the simulation
	@param box_dim (double): Dimensions of the simulation's 3D box
	@param timestep (double): Timestep of the simulation
	@param initial_velocity (double): Max initial velocity of the particles per timestep
	@param sigma (double): Size of the particles in distance units
	@param mass (double): Mass of the particles
*/
class SimulationEngine {
private:
    // Potential
    PotentialField* potential{};

    // Simulation parameters
    int n_particles;
    double box_x, box_y, box_z;
    double timestep;
    int n_steps;
    int current_step = 0;

    // Particle parameters
    double initial_velocity;
    double sigma;
    double mass;



    // Particle Arrays
    particle* particles_t0{};
    particle* particles_t1{};
    particle* particles_t2{};

    // Output
    bool write_output;
    SimulationOutput* output{};

    // Initialize Particles
    void particle_init();

public:

    SimulationEngine(PotentialField* potential, int n_particles, double box_dim, double timestep, int n_step, double initial_velocity, double sigma, double mass);

	// Initialize SimulationEngine Instance
	SimulationEngine(PotentialField* potential, int n_particles, double box_dim, double timestep, int n_step, double initial_velocity, double sigma, double mass, bool write_output);

	// Destructor
	~SimulationEngine();

	int get_n_particles() {
		return this->n_particles;
	}

	int get_current_step() {
		return this->current_step;
	}

	double* get_box_dim();

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
	double nearest_image(particle p1, particle p2);

	/**
	 * @brief Runs the simulation.
	 *
	 * This function runs the simulation for the number of steps specified in the constructor.
	 * It uses the Verlet algorithm to integrate the equations of motion.
	 */
	void run();

};
