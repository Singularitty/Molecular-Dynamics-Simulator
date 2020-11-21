# Molecular-Dynamics-Simulator
Molecular dynamics simulator implemented using the Verlet Algorithm.

Currently has the Lennard Jonnes potencial and a elastic monopole pair potential (equivalent to eletrostatic monopole) implemented. These can be switched by uncomenting sections regarding these potentials in the main cpp file.
This program outputs the following files:
- Outputs particle positions to a positions.xyz file, which can be viewed using a software like OVITO.
- Outputs potencial, kinetic and total energy per timestep in a diferent file called energy.txt.

Example of Lennard-Jones Gas simulation:
- Simulated 108 particles;
- xyz file output viewed in ovito.
![](lennard_jonnes_example.gif)


To do:
- Implement particle charge (multipole momentum for elastic multipole potential);
- Implement particle rotation (important for potentials that depend on angle between two particle directions);
- Implement a easy way to switch potentials;
- Implement a thermostat;
- Parallelize the code to improve performance;

