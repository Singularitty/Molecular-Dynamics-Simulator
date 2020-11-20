//Elastic Multipole Pair Interaction for Nematic Colloids

#include <iostream>
#include <math.h>

using namespace std;

double pi = M_PI;

// Elastic Monopole same momentum (Equivalent to Eletrostatic Monopole)
double U(double r)
{
    double U;
    
    U = (4*pow(b,2)*k*pi*pow(r_eff,2))/r + (4*pow(b,2)*k*pi*(r - rc)*pow(r_eff,2))/pow(rc,2) - (4*pow(b,2)*k*pi*pow(r_eff,2))/rc;

    return U;
}

double f(double pos1[3], double pos2[3], int direction)
{
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


    f = (4*pow(b,2)*k*pi*pow(r_eff,2)*dir)/pow(pow(x,2) + pow(y,2) + pow(z,2),1.5) - 
   (4*pow(b,2)*k*pi*pow(r_eff,2)*dir)/(pow(rc,2)*sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
    
    return f;
}
