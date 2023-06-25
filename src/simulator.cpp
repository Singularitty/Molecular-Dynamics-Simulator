#include <iostream>
#include <string>
#include "Engine/SimulationEngine.hpp"
#include "Parser/SettingsParser.hpp"


int main(int argc, char** argv) {

	if (argc != 2) {
		std::cout << "Usage: ./simulator <settings file>" << std::endl;
		return 1;
	}

    try {
        auto* parser = new SettingsParser(argv[1]);
        parser->parse_sim_params();
        parser->parse_potential_params();
        SimulationEngine* sim = parser->getSimulation();


        std::cout << "Starting simulation..." << std::endl;

        sim->run();

        std::cout << "Simulation finished." << std::endl;

        delete parser;
        delete sim;

        return 0;
    } catch (const std::runtime_error& e) {
        std::cout << "Caught runtime error: " << e.what() << std::endl;
        return 2;
    } catch (...) {
        std::cout << "Caught unknown error" << std::endl;
        return 3;
    }

}