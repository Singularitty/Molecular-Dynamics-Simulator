#pragma once
#include "PotentialField.hpp"

class ElasticMonopole : public PotentialField {
public:
	double k, b, r_eff;

	ElasticMonopole(double box_x, double box_y, double box_z, double k, double b, double r_eff, double cutoff_radius)
		: PotentialField(box_x, box_y, box_z, cutoff_radius), k(k), b(b), r_eff(r_eff) {}

	double potential(double particle_distance);

	double force(particle p1, particle p2, int direction);
};