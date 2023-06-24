// Lennard Jones Potential

#pragma once
#include <iostream>
#include <cmath>

using namespace std;

// Potencial Energy
double U(double ri, double b, double k, double r_eff, double rc);

// Force between two particles
double f(double pos1[3], double pos2[3], int direction, double *box, double b, double k, double r_eff, double rc);
