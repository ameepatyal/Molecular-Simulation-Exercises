# NPT Ensemble Metropolis Monte Carlo

This sub-directory contains a modular C++ implementation of a Metropolis Monte Carlo simulation operating in the isobaric-isothermal (NPT) ensemble.

Unlike the NVT ensemble, this workflow introduces volume fluctuation moves to maintain a constant specified pressure. This serves as a fundamental demonstration of density equilibration algorithms used prior to production-level molecular dynamics.

## System Parameters

*   **Ensemble**: NPT (Constant Number of particles, Pressure, and Temperature)
*   **Particle Count**: 4374 atoms
*   **Initial Simulation Box**: Cubic with an edge length of 18.0 units
*   **Target Reduced Pressure (Pr)**: 0.989
*   **Potential**: Truncated and Shifted Lennard-Jones (TSLJ) with a cutoff distance of 2.5σ
*   **Simulation Length**: 2,000,000 Monte Carlo steps

## Code Structure

The codebase is refactored to adhere to the DRY (Don't Repeat Yourself) principle, separating core physics calculations from execution scripts.

*   `include/physics_utils.h`: Header file containing declarations for periodic boundary condition (PBC) distance calculations and the TSLJ potential.
*   `src/physics_utils.cpp`: The implementation of the shared physics calculations.
*   `src/03_npt_lattice_setup.cpp`: Generates 16x17x17 initial grid, establishes starting coordinates, and computes initial Radial Distribution Function (RDF).
*   `src/04_mc_npt_simulation.cpp`: Reads the initialized coordinates and executes Monte Carlo sampling loop, balancing translational particle displacements with simulation box volume scaling.

## Compilation & Execution

Compile the setup script and run it to generate the initial lattice and starting RDF:

```bash
g++ -I./include src/03_npt_lattice_setup.cpp src/physics_utils.cpp -o setup_npt
./setup_npt