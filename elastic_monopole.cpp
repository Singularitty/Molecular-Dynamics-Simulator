#include "elastic_monopole.h"

double pi = 3.1415926535;

double U(double r, double b, double k, double r_eff, double rc) {
    double U = (4*pow(b, 2) * k * pi * pow(r_eff, 2))/r + (4 * pow(b, 2) * k * pi * (r - rc) *
            pow(r_eff, 2))/pow(rc, 2) - (4 * pow(b, 2) * k * pi * pow(r_eff, 2)) / rc;
    return U;
}

double f(double pos1[3], double pos2[3], int direction, double *box, double b, double k, double r_eff, double rc) {
    double x, y, z, dir = 0;
    x = pos1[0] - pos2[0];
    x = x - box[0] * round(x / box[0]);
    y = pos1[1] - pos2[1];
    y = y - box[1] * round(y / box[1]);
    z = pos1[2] - pos2[2];
    z = z - box[2] * round(z/box[2]);
    if(direction == 0)
        dir = x;
    if(direction == 1)
        dir = y;
    if(direction == 2)
        dir = z;

    double f = (4 * pow(b, 2) * k * pi * pow(r_eff, 2) * dir) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5) -
        (4 * pow(b, 2) * k * pi * pow(r_eff, 2) * dir) / (pow(rc, 2) * sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)));
    return f;
}
