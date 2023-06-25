#include "ElasticMonopole.hpp"

double ElasticMonopole::potential(double particle_distance) {

    double u;

    u = (4 * pow(b, 2) * k * pi * pow(r_eff, 2)) / particle_distance + (4 * pow(b, 2) * k * pi * (particle_distance - cutoff_radius) * pow(r_eff, 2)) / pow(cutoff_radius, 2) - (4 * pow(b, 2) * k * pi * pow(r_eff, 2)) / cutoff_radius;

    return u;
}

double ElasticMonopole::force(particle p1, particle p2, int direction) {

    double f;
    double x, y, z;
    double dir = 0;

    x = p1.x - p2.x;
    x = x - box_x * round(x / box_x);
    y = p1.y - p2.y;
    y = y - box_y * round(y / box_y);
    z = p1.z - p2.z;
    z = z - box_z * round(z / box_z);

    if (direction == 0) {
        dir = x;
    }
    else if (direction == 1) {
        dir = y;
    }
    else if (direction == 2) {
        dir = z;
    }

    f = (4 * pow(b, 2) * k * pi * pow(r_eff, 2) * dir) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5) -
            (4 * pow(b, 2) * k * pi * pow(r_eff, 2) * dir) / (pow(cutoff_radius, 2) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));

    return f;
}
