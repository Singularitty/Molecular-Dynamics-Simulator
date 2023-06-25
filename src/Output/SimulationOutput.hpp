#pragma once
#include <fstream>
#include <string>
#include <chrono>
#include <sstream>
#include <iomanip>
#include "../Engine/particle.hpp"

class SimulationEngine;

class SimulationOutput {
private:
	SimulationEngine* sim;
	std::unique_ptr<std::ofstream> positions;
	std::unique_ptr<std::ofstream> energy;

	// Persistent SimulationEngine Data
	double box_x, box_y, box_z;
	int n_particles;
public:

	explicit SimulationOutput(SimulationEngine* sim);

    ~SimulationOutput() {
		if (positions->is_open()) {
			positions->close();
		}
		if (energy->is_open()) {
			energy->close();
		}
	}

	void write_energy(double total_energy, double potential_energy, double kinectic_energy);

	void write_positions_frame(particle* particles);
};