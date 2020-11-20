// Lennard Jones Potential

#include <iostream
#include <math.h>

using namespace std;

// Potencial Energy
double U(double r) {

    double U;
    U = 4.*epsilon*(-(pow(sigma,6)/pow(r,6)) + pow(sigma,12)/pow(r,12)) - 4.*epsilon*(r - rc)*((6.*pow(sigma,6))/pow(rc,7) - (12.*pow(sigma,12))/pow(rc,13)) - 
    4.*epsilon*(-(pow(sigma,6)/pow(rc,6)) + pow(sigma,12)/pow(rc,12));
    
    return U;
}

// Force between two particles
double f(double pos1[3], double pos2[3], int direction) {

    double f;
    double x,y,z, dir;
    x = pos1[0] - pos2[0];
    x = x - box[0]*round(x/box[0]);
    y = pos1[1] - pos2[1];
    y = y - box[1]*round(y/box[1]);
    z = pos1[2] - pos2[2];
    z = z - box[2]*round(z/box[2]);
    if (direction == 0) dir = x;
    if (direction == 1) dir = y;
    if (direction == 2) dir = z;

    f = -1.*((-4.*epsilon*((6.*pow(sigma,6))/pow(rc,7) - (12.*pow(sigma,12))/pow(rc,13))*dir)/sqrt(pow(x,2) + pow(y,2) + pow(z,2)) + 
    4.*epsilon*((-12.*pow(sigma,12)*dir)/pow(pow(x,2) + pow(y,2) + pow(z,2),7) + (6.*pow(sigma,6)*dir)/pow(pow(x,2) + pow(y,2) + pow(z,2),4)));

    return f;
}