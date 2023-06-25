#include "SettingsParser.hpp"

SettingsParser::SettingsParser(const std::string& settings_file) {
    this->settings_file = std::make_unique<std::ifstream>(settings_file);


    if (!this->settings_file->is_open()) {
        throw std::runtime_error("Could not open settings file");
    }
}


void SettingsParser::parse_sim_params() {

    int found_settings = 0;

    std::string line;
    while (std::getline(*settings_file, line)) {
        std::istringstream is_line(line);
        std::string key;
        if (std::getline(is_line, key, '=')) {
            std::string value;
            if (std::getline(is_line, value)) {
                //std::cout << "Key: " << key << " Value: " << value << std::endl;
                if (key == "PARTICLES") {
                    this->particles = std::stoi(value);
                    found_settings++;
                }
                else if (key == "BOX") {
                    this->box = std::stod(value);
                    found_settings++;
                }
                else if (key == "STEPS") {
                    this->steps = std::stoi(value);
                    found_settings++;
                }
                else if (key == "TIMESTEP") {
                    this->timestep = std::stod(value);
                    found_settings++;
                }
                else if (key == "DISPLACEMENT") {
                    this->displacement = std::stod(value);
                    found_settings++;
                }
                else if (key == "PARTICLE_SIZE") {
                    this->particle_size = std::stod(value);
                    found_settings++;
                }
                else if (key == "MASS") {
                    this->mass = std::stod(value);
                    found_settings++;
                }
            }
        }
    }

    if (found_settings != 7) {
        throw std::runtime_error("Not all simulation parameters were found");
    }

    std::cout << "Simulation parameters:" << std::endl;
    std::cout << " - Particles: " << this->particles << std::endl;
    std::cout << " - Box: " << this->box << std::endl;
    std::cout << " - Steps: " << this->steps << std::endl;
    std::cout << " - Timestep: " << this->timestep << std::endl;
    std::cout << " - Displacement: " << this->displacement << std::endl;
    std::cout << " - Particle size: " << this->particle_size << std::endl;
    std::cout << " - Mass: " << this->mass << std::endl;

}


void SettingsParser::parse_potential_params() {

    // Set file pointer to the beginning of the file
    settings_file->clear();
    settings_file->seekg(0, std::ios::beg);

    // Determine potential type
    std::string potential_type;
    std::string line;
    while (std::getline(*settings_file, line)) {
        std::istringstream is_line(line);
        std::string key;
        if (std::getline(is_line, key, '=')) {
            std::string value;
            if (std::getline(is_line, value)) {
                if (key == "POTENTIAL") {
                    potential_type = value;
                    break;
                }
            }
        }
    }

    // Set file pointer to the beginning of the file
    settings_file->clear();
    settings_file->seekg(0, std::ios::beg);

    // Remove pesky newline character
    // This is bad practice, but it works
    potential_type.pop_back();

    // Parse potential specific parameters
    if (potential_type == "LennardJones") {
        double dispersion, cutoff;

        while (std::getline(*settings_file, line)) {
            std::istringstream is_line(line);
            std::string key;
            if (std::getline(is_line, key, '=')) {
                std::string value;
                if (std::getline(is_line, value)) {
                    if (key == "DISPERSION_ENERGY") {
                        dispersion = std::stod(value);
                    }
                    else if (key == "CUTOFF") {
                        cutoff = std::stod(value);
                    }
                }
            }
        }

        potential = new LennardJones(this->box, this->box, this->box, dispersion, this->particle_size, cutoff);
        sim = new SimulationEngine(potential, this->particles, this->box, this->timestep, this->steps, this->displacement, this->particle_size, this->mass);

        std::cout << "Lennard Jones Potential parameters:" << std::endl;
        std::cout << " - Dispersion: " << dispersion << std::endl;
        std::cout << " - Cutoff: " << cutoff << std::endl;

    }
    else if (potential_type == "ElasticMonopole") {
        double momentum, effective_radius, cutoff, k;

        while (std::getline(*settings_file, line)) {
            std::istringstream is_line(line);
            std::string key;
            if (std::getline(is_line, key, '=')) {
                std::string value;
                if (std::getline(is_line, value)) {
                    if (key == "MOMENTUM") {
                        momentum = std::stod(value);
                    }
                    else if (key == "EFFECTIVE_RADIUS") {
                        effective_radius = std::stod(value);
                    }
                    else if (key == "CUTOFF") {
                        cutoff = std::stod(value);
                    }
                    else if (key == "K") {
                        k = std::stod(value);
                    }
                }
            }
        }

        potential = new ElasticMonopole(this->box, this->box, this->box, k, momentum, effective_radius, cutoff);
        sim = new SimulationEngine(potential, this->particles, this->box, this->timestep, this->steps, this->displacement, this->particle_size, this->mass);

        std::cout << "Elastic Monopole Potential parameters:" << std::endl;
        std::cout << " - Momentum: " << momentum << std::endl;
        std::cout << " - Effective radius: " << effective_radius << std::endl;
        std::cout << " - Cutoff: " << cutoff << std::endl;
        std::cout << " - K: " << k << std::endl;


    }
    else {
        throw std::runtime_error("Unknown Potential");
    }

}
