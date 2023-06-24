#pragma once
#include <iostream>
#include <cmath>

using namespace std;

// Elastic Monopole same momentum (Equivalent to Eletrostatic Monopole)
double U(double r, double b, double k , double r_eff, double rc);

double f(double pos1[3], double pos2[3], int direction, double *box, double b, double k, double r_eff, double rc);
