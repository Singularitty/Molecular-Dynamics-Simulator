#include "SimulationOutput.hpp"
#include "../Engine/SimulationEngine.hpp"

SimulationOutput::SimulationOutput(SimulationEngine* sim)
        : sim(sim)
{

    // Get Timestamp
    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);

    std::stringstream timestamp;
    timestamp << std::put_time(std::localtime(&now_c), "%Y-%m-%d-%H-%M-%S");

    // Create filenames
    std::string positions_filename = "positions-" + timestamp.str() + ".xyz";
    std::string energy_filename = "energy-" + timestamp.str() + ".csv";

    // Open Position xyz file

    positions = std::make_unique<std::ofstream>(positions_filename);

    if (!positions->is_open()) {
        throw std::runtime_error("Could not open file " + positions_filename);
    }

    // Open Energy csv file

    energy = std::make_unique<std::ofstream>(energy_filename);

    if (!energy->is_open()) {
        throw std::runtime_error("Could not open file " + energy_filename);
    }

    // Write header to energy file
    *energy << "Step,Total Energy,Potential Energy,Kinetic Energy" << std::endl;


    // Get Constant Simulation Parameters
    n_particles = sim->get_n_particles();

    double* box_dims = sim->get_box_dim();
    box_x = box_dims[0];
    box_y = box_dims[1];
    box_z = box_dims[2];
    delete[] box_dims;
}

void SimulationOutput::write_energy(double total_energy, double potential_energy, double kinectic_energy) {
    *energy << sim->get_current_step() << "," << total_energy << ","  << potential_energy << "," << kinectic_energy << std::endl;
}

void SimulationOutput::write_positions_frame(particle* particles) {

    *positions << "ITEM: TIMESTEP" << std::endl;
    *positions << sim->get_current_step() << std::endl;
    *positions << "ITEM: NUMBER OF ATOMS" << std::endl;
    *positions << n_particles << std::endl;
    *positions << "ITEM: BOX BOUNDS pp pp pp" << std::endl;
    *positions << 0 << "\t" << box_x << std::endl;
    *positions << 0 << "\t" << box_y << std::endl;
    *positions << 0 << "\t" << box_z << std::endl;
    *positions << "ITEM: ATOMS id x y z" << std::endl;
    for (int i = 0; i < this->n_particles; i++) {
        *positions << i << "\t" << particles[i].x << "\t" << particles[i].y << "\t" << particles[i].z << std::endl;
    }
}
