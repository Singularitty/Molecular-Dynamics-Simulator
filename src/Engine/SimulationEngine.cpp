#include "SimulationEngine.hpp"

SimulationEngine::SimulationEngine(PotentialField* potential, int n_particles, double box_dim, double timestep, int n_step, double initial_velocity, double sigma, double mass)
    : SimulationEngine(potential, n_particles, box_dim, timestep, n_step, initial_velocity, sigma, mass, true) {}

// Initialize Simulation Instance
SimulationEngine::SimulationEngine(PotentialField* potential, int n_particles, double box_dim, double timestep, int n_step, double initial_velocity, double sigma, double mass, bool write_output)
        : potential(potential), n_particles(n_particles), box_x(box_dim), box_y(box_dim), box_z(box_dim), n_steps(n_step), initial_velocity(initial_velocity), sigma(sigma), mass(mass), write_output(write_output)
{
    // Particle Array Initialization for 3 Time Instances

    // Time Instance 0
    particles_t0 = new particle[this->n_particles];
    // Time Instance 1
    particles_t1 = new particle[this->n_particles];
    // Time Instance 2
    particles_t2 = new particle[this->n_particles];

    // Create Output
    if (write_output) {
        output = new SimulationOutput(this);
    }

}

// Initialize Particles
void SimulationEngine::particle_init() {

    bool colision;

    // Random Number from Hardware
    std::random_device dev;
    // Seed Pesudo-Random Generator
    std::mt19937 gen(dev());
    // Real distribution from 0 to 1
    std::uniform_real_distribution<> distr(0, 1);

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

        // While there is a colision with another particle generate a new one
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

// Destructor
SimulationEngine::~SimulationEngine() {
    delete[] particles_t0;
    delete[] particles_t1;
    delete[] particles_t2;

    if (write_output) {
        delete output;
    }
}


double* SimulationEngine::get_box_dim() {
    double* box_dim = new double[3];
    box_dim[0] = this->box_x;
    box_dim[1] = this->box_y;
    box_dim[2] = this->box_z;
    return box_dim;
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
double SimulationEngine::nearest_image(particle p1, particle p2) {

    double x, y, z;

    x = p1.x - p2.x;
    x = x - this->box_x * round(x / this->box_x);
    y = p1.y - p2.y;
    y = y - this->box_y * round(y / this->box_y);
    z = p1.z - p2.z;
    z = z - this->box_z * round(z / this->box_z);

    return sqrt(x * x + y * y + z * z);
}

/**
 * @brief Runs the simulation.
 *
 * This function runs the simulation for the number of steps specified in the constructor.
 * It uses the Verlet algorithm to integrate the equations of motion.
 */
void SimulationEngine::run() {

    double particle_distance = 0.0;
    double p_energy, k_energy, t_energy;
    double vel, acc;

    std::cout << "Generating Particles..." << std::endl;

    // Initialize Random Particle Positions
    this->particle_init();

    std::cout << "Running Simulation..." << std::endl;

    for (int step = 0; step < this->n_steps; step++) {

        this->current_step = step;

        p_energy = 0.;
        k_energy = 0.;
        t_energy = 0.;

        for (int i = 0; i < this->n_particles; i++) {
            vel = 0.;
            // For each space dimensions
            for (int dim = 0; dim < 3; dim++) {

                acc = 0.;
                for (int j = 0; j < this->n_particles; j++) {

                    particle_distance = this->nearest_image(particles_t1[i], particles_t1[j]);

                    if (i != j && particle_distance < potential->cutoff_radius) {
                        acc += potential->force(particles_t0[i], particles_t0[j], dim);
                    }
                }

                switch (dim)
                {
                    case 0:
                        this->particles_t2[i].x = fmod(2. * particles_t1[i].x - particles_t0[i].x + (acc / this->mass) * this->timestep * this->timestep + 2. * this->box_x, this->box_x);
                        break;
                    case 1:
                        this->particles_t2[i].y = fmod(2. * particles_t1[i].y - particles_t0[i].y + (acc / this->mass) * this->timestep * this->timestep + 2. * this->box_y, this->box_y);
                        break;
                    case 2:
                        this->particles_t2[i].z = fmod(2. * this->particles_t1[i].z - this->particles_t0[i].z + (acc / this->mass) * this->timestep * this->timestep + 2. * this->box_z, this->box_z);
                        break;
                }
            }

            for (int j = i + 1; j < this->n_particles; j++) {

                particle_distance = this->nearest_image(particles_t1[i], particles_t1[j]);
                if (particle_distance < potential->cutoff_radius) {
                    p_energy += potential->potential(particle_distance);
                }
            }

            vel = this->nearest_image(particles_t2[i], particles_t0[i]) / (2. * this->timestep);
            k_energy += 0.5 * this->mass * vel * vel;

            particles_t0[i] = particles_t1[i];
            particles_t1[i] = particles_t2[i];

        }

        t_energy = p_energy + k_energy;

        if (write_output) {
            output->write_energy(t_energy, p_energy, k_energy);
            output->write_positions_frame(this->particles_t0);
        }
    }

}
