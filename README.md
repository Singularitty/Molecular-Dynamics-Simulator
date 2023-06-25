# Molecular-Dynamics-Simulator
Molecular dynamics simulator implemented using the Verlet Algorithm.

Currently has the Lennard Jones potential and a elastic monopole pair potential (equivalent to electrostatic monopole) implemented.
This program outputs the following files for each run:
- Outputs particle positions to a positions.xyz file, which can be viewed using a software like OVITO.
- Outputs potential, kinetic and total energy per timestep in a different file called energy.txt.

Output files are timestamped with datetime for easier identification between runs.

## Usage

To compile the program run make in the repository's directory
```
make
```

This will generate a simulator binary file which can then ran using

```
./simulator <settings file>
```

The settings file contains the simulation parameters and potential to be used alongside it's parameters. Check the example settings files in the repository for the Lennard-Jones and Elastic Monopole potentials.

## Example of Lennard-Jones Gas simulation:
- Simulated 108 particles;
- xyz file output viewed in ovito.
![](lennard_jonnes_example.gif)



#### I won't be improving this simulator any further (for now at least), but here are some features that could be implemented:

- More potentials
- Particle charge (multipole momentum for elastic multipole potential);
- Particle rotation (important for potentials that depend on angle between two particle directions);
- A thermostat;
- Parallelization of the code to improve performance;

