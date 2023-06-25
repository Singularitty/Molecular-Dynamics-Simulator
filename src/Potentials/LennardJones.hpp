#pragma once
#include "PotentialField.hpp"

class LennardJones : public PotentialField {
public:

	double epsilon, sigma;

	LennardJones(double box_x, double box_y, double box_z, double epsilon, double sigma, double cutoff_radius)
		: PotentialField(box_x, box_y, box_z, cutoff_radius), epsilon(epsilon), sigma(sigma) {}

	double potential(double particle_distance);

	double force(particle p1, particle p2, int direction);
};
