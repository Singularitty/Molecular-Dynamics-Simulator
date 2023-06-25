#pragma once
#include <cmath>

// This Struct represents a point particle in 3D space
struct particle {
    double x;
    double y;
    double z;

    particle() : x(0), y(0), z(0) {}

    particle(double x, double y, double z)
            : x(x), y(y), z(z) {};
};

// Distance between two particles
double distance(particle p1, particle p2);