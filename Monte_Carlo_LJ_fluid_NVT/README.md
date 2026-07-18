# NVT Ensemble Metropolis Monte Carlo

This sub-directory contains a modular C++ implementation of a Metropolis Monte Carlo simulation operating in the canonical (NVT) ensemble. 

The workflow establishes an initial simple cubic lattice and samples configuration space to compute thermodynamic properties and structural features while utilizing a Truncated and Shifted Lennard-Jones (TSLJ) potential.

## System Parameters
*   **Ensemble**: NVT (Constant Number of particles, Volume, and Temperature)
*   **Particle Count**: 2700 atoms
*   **Simulation Box**: Cubic, with an edge length of 15.0 units
*   **Potential**: Lennard-Jones with a cutoff distance of 2.5 $\sigma$
*   **Simulation Length**: 10,000 Monte Carlo steps

## Code Structure
The codebase is refactored to adhere to the DRY (Don't Repeat Yourself) principle, separating core physics calculations from execution scripts.

*   `include/physics_utils.h`: Header file containing the declarations for periodic boundary condition (PBC) distance calculations, TSLJ potential, and total energy calculations.
*   `src/physics_utils.cpp`: The implementation of the shared physics calculations.
*   `src/01_lattice_setup.cpp`: Generates the 14x14x14 grid, removes coordinates to reach the target number density, and outputs the initial configuration.
*   `src/02_mc_nvt_simulation.cpp`: Reads initial coordinates, executes Metropolis sampling algorithm, and evaluates radial distribution function (RDF).

## Compilation & Execution

Compile the setup script and run it to generate the initial lattice:
```bash
g++ -I./include src/01_lattice_setup.cpp src/physics_utils.cpp -o setup_lattice
./setup_lattice
