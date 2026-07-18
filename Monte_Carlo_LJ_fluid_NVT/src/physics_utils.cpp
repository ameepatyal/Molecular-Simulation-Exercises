#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>

// Function to compute the distance between two points with periodic boundary conditions
double PBC_Distance(double x1, double y1, double z1, double x2, double y2, double z2, double boxSize) {
    double dx = std::abs(x1 - x2);
    double dy = std::abs(y1 - y2);
    double dz = std::abs(z1 - z2);

    // Apply periodic boundary conditions
    if (dx > boxSize / 2.0) dx = boxSize - dx;
    if (dy > boxSize / 2.0) dy = boxSize - dy;
    if (dz > boxSize / 2.0) dz = boxSize - dz;

    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Function to compute the truncated and shifted Lennard-Jones potential with a cutoff of 2.5 sigma
double LJ_Potential(double r, double epsilon, double sigma, double cutoff) {
    if (r > cutoff) {
        return 0.0;
    } else {
        double sigmaOverR = sigma / r;
        double sigmaOverR6 = sigmaOverR * sigmaOverR * sigmaOverR * sigmaOverR * sigmaOverR * sigmaOverR;
        double sigmaOverR6cutoff = (sigma / cutoff) * (sigma / cutoff)* (sigma / cutoff)* (sigma / cutoff)* (sigma / cutoff)* (sigma / cutoff);
        double sigmaOverR12 = sigmaOverR6 * sigmaOverR6;
        double sigmaOverR12cutoff = sigmaOverR6cutoff*sigmaOverR6cutoff;
        double potential = 4.0 * epsilon * (sigmaOverR12 - sigmaOverR6);
        double potential_cutoff = 4.0 * epsilon *( sigmaOverR12cutoff - sigmaOverR6cutoff );
        return potential - potential_cutoff;
    }
}

// Function to calculate the total energy of the system
double CalculateTotalEnergy(const std::vector<double>& xCoordinates, const std::vector<double>& yCoordinates,
                            const std::vector<double>& zCoordinates, double epsilon, double sigma, double cutoff, double boxSize) {
    double totalEnergy = 0.0;
    int numParticles = xCoordinates.size();

    for (int i = 0; i < numParticles; ++i) {
        for (int j = i + 1; j < numParticles; ++j) {
            double distance = PBC_Distance(xCoordinates[i], yCoordinates[i], zCoordinates[i],
                                              xCoordinates[j], yCoordinates[j], zCoordinates[j], boxSize);
            double potential = LJ_Potential(distance, epsilon, sigma, cutoff);
            totalEnergy += potential;
        }
    }

    return totalEnergy;
}