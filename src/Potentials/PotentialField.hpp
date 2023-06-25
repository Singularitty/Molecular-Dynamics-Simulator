#pragma once
#include <cmath>
#include "../Engine/particle.hpp"

const double pi = 3.14159265358979323846;

class PotentialField {
public:
    double box_x, box_y, box_z;
    double cutoff_radius;

    PotentialField(double box_x, double box_y, double box_z, double cutoff_radius)
            : box_x(box_x), box_y(box_y), box_z(box_z), cutoff_radius(cutoff_radius) {}

    virtual double force(particle p1, particle p2, int dir) = 0;

    virtual double potential(double particle_distance) = 0;
};