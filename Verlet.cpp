#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>

using namespace std;

// Simulation Parameters
int N = 4*pow(3,3);       // Number of particles
double box[3] = {10.0,10.0,10.0};  // Dimensions of Simulation Box
double timestep = 1.e-5 / 20.;     // Timestep
int Num_Steps = 1e+5;              // Number of Timesteps
double inicial_max_displacement = 1.e-5; // Determines the inicial velocity of the particles
double sigma = 1.0;         // Size of the particle (distance units)
double m = 1.0;             // Mass

/* 
// Lennard Jonnes Potential
double epsilon = 1.0;       // Dispersion energy (energy units)
double rc = 2.5;            // Cutoff radius (distance units)
#include "Lennard_Jonnes.h"
 */


// Elastic Multipole 
double b = 1.;    // Momentum
double r_eff = 1.;    // Efective radius
double rc = 3.;     // Cutoff radius
double k = 1.;
#include "elastic_monopole.h"

// Nearest Image Convention
double NIC(double pos1[3],double pos2[3]) {
    
    double x, y, z;

    x = pos1[0] - pos2[0];
    x = x - box[0]*round(x/box[0]);
    y = pos1[1] - pos2[1];
    y = y - box[1]*round(y/box[1]);
    z = pos1[2] - pos2[2];
    z = z - box[2]*round(z/box[2]);

    return sqrt(x*x + y*y + z*z);
}

// Distance Between two points
double r(double pos1[3], double pos2[3]) {
    
    double r;
    r = sqrt(pow((pos1[0]-pos2[0]),2)+pow((pos1[1]-pos2[1]),2)+pow((pos1[2]-pos2[2]),2));

    return r;
}

int main() {

    // List of particle coordinates in 3D space in 3 diferent time instances
    double ***Pos = new double**[3];         // Instance in time (0 = t-delta_t, 1 = t, 2 = t + delta_t)
    for (int i = 0; i < 3; i++)
    {
        Pos[i] = new double*[N];            // Particle id

        for (int j = 0; j < N; j++)
            Pos[i][j] = new double[3];      // Coordinate in space (0 = x, 1 = y, 2 = z)
    }

    // Particle Generator
    bool Colision;

    cout << "Generating Particles...\n";

    for (int i = 0; i < N; i++)
    {
        if (i == 0)
        {
            Pos[0][i][0] = drand48()*box[0];
            Pos[0][i][1] = drand48()*box[1];
            Pos[0][i][2] = drand48()*box[2];

            Pos[1][i][0] = inicial_max_displacement*(1.-2.*drand48()) + Pos[0][i][0];
            Pos[1][i][1] = inicial_max_displacement*(1.-2.*drand48()) + Pos[0][i][1];
            Pos[1][i][2] = inicial_max_displacement*(1.-2.*drand48()) + Pos[0][i][2];
        }
        else
        {
            Colision = true;
            while (Colision == true)
            {
                Colision = false;

                Pos[0][i][0] = drand48()*box[0];
                Pos[0][i][1] = drand48()*box[1];
                Pos[0][i][2] = drand48()*box[2];
                
                for (int j = 0; j < i; j++)
                {
                    if (NIC(Pos[0][i],Pos[0][j]) < 1.5*sigma)
                        Colision = true;
                }
            }
    
            Pos[1][i][0] = inicial_max_displacement*(1.-2.*drand48()) + Pos[0][i][0];
            Pos[1][i][1] = inicial_max_displacement*(1.-2.*drand48()) + Pos[0][i][1];
            Pos[1][i][2] = inicial_max_displacement*(1.-2.*drand48()) + Pos[0][i][2];
        }
    }

    cout << "Particle Generation Complete. \n";

    ofstream positions;
    ofstream energy;

    positions.open("positions.xyz");
    energy.open("energy.txt");

    positions << "ITEM: TIMESTEP\n" << 0 << endl;
    positions << "ITEM: NUMBER OF ATOMS\n" << N << endl;
    positions << "ITEM: BOX BOUNDS pp pp pp" << endl;
    positions << 0 << "\t" << box[0] << endl;
    positions << 0 << "\t" << box[1] << endl;
    positions << 0 << "\t" << box[2] << endl;
    positions << "ITEM: ATOMS id x y z" << endl;
    for (int i = 0; i < N; i++)
        {
            positions << i << "\t" << Pos[0][i][0] << "\t" << Pos[0][i][1] << "\t" << Pos[0][i][2] << endl;
        }
    positions << "ITEM: TIMESTEP\n" << 1 << endl;
    positions << "ITEM: NUMBER OF ATOMS\n" << N << endl;
    positions << "ITEM: BOX BOUNDS pp pp pp" << endl;
    positions << 0 << "\t" << box[0] << endl;
    positions << 0 << "\t" << box[1] << endl;
    positions << 0 << "\t" << box[2] << endl;
    positions << "ITEM: ATOMS id x y z" << endl;
    for (int i = 0; i < N; i++)
        {
            positions << i << "\t" << Pos[1][i][0] << "\t" << Pos[1][i][1] << "\t" << Pos[1][i][2] << endl;
        }

    // Verlet Algorithm

    double a, v, Potencial_energy, Kinetic_energy, Total_energy;
    double particle_distance;

    for (int t = 0; t < Num_Steps; t++) {

        Potencial_energy = 0.;
        Kinetic_energy = 0.;
        Total_energy = 0.;

        for (int i = 0; i < N; i++) 
        {
            v = 0.;

            for (int k = 0; k < 3; k++)
            {
                a = 0.;

                for (int j = 0; j < N; j++)
                {
                    particle_distance = NIC(Pos[1][i],Pos[1][j]);
                    if (j != i && particle_distance < rc)
                    {
                        a += f(Pos[1][i],Pos[1][j],k);
                    }
                }

                Pos[2][i][k] = fmod(2.*Pos[1][i][k] - Pos[0][i][k] + (a/m)*timestep*timestep + 2.*box[k],box[k]);

            }

            for (int j = i+1; j < N; j++)
            {
                particle_distance = NIC(Pos[1][i],Pos[1][j]);
                if (particle_distance < rc)
                {
                    Potencial_energy += U(particle_distance);
                }
            }

            v = NIC(Pos[2][i],Pos[0][i]) / (2.*timestep);
            Kinetic_energy += 0.5*m*v*v;


            for (int k = 0; k < 3; k++)
            {
                Pos[0][i][k] = Pos[1][i][k];
                Pos[1][i][k] = Pos[2][i][k];
            }

        }

        Total_energy = Potencial_energy + Kinetic_energy;
        if (t % 100 == 0)
            energy << t+1 << "\t" << Potencial_energy << "\t" << Kinetic_energy << "\t" << Total_energy << endl;
        
        positions << "ITEM: TIMESTEP\n" << t+2 << endl;
        positions << "ITEM: NUMBER OF ATOMS\n" << N << endl;
        positions << "ITEM: BOX BOUNDS pp pp pp" << endl;
        positions << 0 << "\t" << box[0] << endl;
        positions << 0 << "\t" << box[1] << endl;
        positions << 0 << "\t" << box[2] << endl;
        positions << "ITEM: ATOMS id x y z" << endl;
        for (int i = 0; i < N; i++)
        {
            positions << i << "\t" << Pos[2][i][0] << "\t" << Pos[2][i][1] << "\t" << Pos[2][i][2] << endl;
        }
    }

    return 0;
}
