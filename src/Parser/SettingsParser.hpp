#pragma once
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "../Engine/SimulationEngine.hpp"
#include "../Potentials/PotentialField.hpp"
#include "../Potentials/LennardJones.hpp"
#include "../Potentials/ElasticMonopole.hpp"

class SettingsParser {
private:
	std::unique_ptr<std::ifstream> settings_file;
	SimulationEngine* sim;
	PotentialField* potential;

	// SimulationEngine parameters
	int particles, steps;
	double box, timestep, displacement, particle_size, mass;
public:

	explicit SettingsParser(const std::string& settings_file) ;

	~SettingsParser() {
		if (settings_file->is_open()) {
			settings_file->close();
		}
	}

	SimulationEngine* getSimulation() {
		return sim;
	}

	void parse_sim_params();

	void parse_potential_params();
	
};