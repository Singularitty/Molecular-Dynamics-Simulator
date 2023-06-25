#include "LennardJones.hpp"

double LennardJones::potential(double particle_distance) {

    double u;

    u = 4.*epsilon*(-(pow(sigma,6)/pow(particle_distance,6)) + pow(sigma,12)/pow(particle_distance,12))
        - 4.*epsilon*(particle_distance - cutoff_radius)*((6.*pow(sigma,6))/pow(cutoff_radius,7) - (12.*pow(sigma,12))/pow(cutoff_radius,13))
        - 4.*epsilon*(-(pow(sigma,6)/pow(cutoff_radius,6)) + pow(sigma,12)/pow(cutoff_radius,12));

    return u;
}

double LennardJones::force(particle p1, particle p2, int direction) {

    double f;
    double x, y, z;
    double dir = 0;

    x = p1.x - p2.x;
    x = x - box_x*round(x/box_x);
    y = p1.y - p2.y;
    y = y - box_y*round(y/box_y);
    z = p1.z - p2.z;
    z = z - box_z*round(z/box_z);

    if (direction == 0) {
        dir = x;
    } else if (direction == 1) {
        dir = y;
    } else if (direction == 2) {
        dir = z;
    }

    f = -1. * ((-4. * epsilon * ((6. * pow(sigma, 6)) / pow(cutoff_radius, 7) - (12. * pow(sigma, 12)) / pow(cutoff_radius, 13)) * dir) / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) +
               4. * epsilon * ((-12. * pow(sigma, 12) * dir) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 7) + (6. * pow(sigma, 6) * dir) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 4)));

    return f;
}

