#pragma once

// Number of particles
// Possible values: x^3 * 4 = 108 256 500 864 1372 2048 2916 4000 5324 6912 8788 10976
int N = 108;

// Dimensions of simulation box
double box[3] = {10.0, 10.0, 10.0};

// Timestep
double timestep = 1.e-5 / 20.0;

// Number of Timesteps
int Num_Steps = 1e+5;

// Determines the initial velocity of the particles
double inicial_max_displacement = 1.e-5;

// Size of the particle (distance units)
double sigma = 1.0;

// Mass
double m = 1.0;

/* Leonard Jonnes Potential
// Dispersion energy (energy units)
double epsilon = 1.0;

// Cutoff radius (distance units)
double rc = 2.5;
*/

/* Elastic Multipole */
// Momentum
double b = 1.0;

// Efective radius
double r_eff = 1.0;

// Cutoff radius
double rc = 3.0;

double k = 1.0;
