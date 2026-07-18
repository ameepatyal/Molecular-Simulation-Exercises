#ifndef PHYSICS_UTILS_H
#define PHYSICS_UTILS_H

#include <vector>

// Computes the distance between two points with periodic boundary conditions
double PBC_Distance(double x1, double y1, double z1, double x2, double y2, double z2, double boxSize);

// Computes the truncated and shifted Lennard-Jones potential
double LJ_Potential(double r, double epsilon, double sigma, double cutoff);

// Calculates the total energy of the system
double CalculateTotalEnergy(const std::vector<double>& xCoordinates, const std::vector<double>& yCoordinates,
                            const std::vector<double>& zCoordinates, double epsilon, double sigma, double cutoff, double boxSize);

#endif