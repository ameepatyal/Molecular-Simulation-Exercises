#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include "physics_utils.h"

int main() {
    double a = 15.0; // Box edge length
    int totalCoordinates = 2744;
    int Gridintervals = 14; // Number of intervals
    double Grid_width = a / Gridintervals;

    std::ofstream outputFile("Grid_coordinates.txt");
    int f = 0;
    if (outputFile.is_open()) {
        for (int i = 0; i < Gridintervals; ++i) {
            for (int j = 0; j < Gridintervals; ++j) {
                for (int k = 0; k < Gridintervals; ++k) {
                    {if (f<2700){
                    double x = i * Grid_width;
                    double y = j * Grid_width;
                    double z = k * Grid_width;
                    outputFile << "H  " << x << " " << y << " " << z << std::endl;
                    ++f;}
                    }
                }
            }
        }

        outputFile.close();
        std::cout << "Grid coordinates have been saved to 'Grid_coordinates.txt'" << std::endl;

    // Define parameters for TSLJ potential calculation
    double boxSize = a;
    double epsilon = 1.0; // LJ potential parameter (adjust as needed)
    double sigma = 1.0;   // LJ potential parameter (adjust as needed)
    double cutoff = 2.5 * sigma; // Cutoff for TSLJ potential
    double temperature = 1.0; // Temperature

    // Read coordinates from the generated file
    std::ifstream inputFile("Coordinates_Restart.txt");

    std::vector<double> xCoordinates;
    std::vector<double> yCoordinates;
    std::vector<double> zCoordinates;

    std::string line;
    while (std::getline(inputFile, line)) {
        double x, y, z;
        if (std::sscanf(line.c_str(), "H %lf %lf %lf", &x, &y, &z) == 3) {
            xCoordinates.push_back(x);
            yCoordinates.push_back(y);
            zCoordinates.push_back(z);
        }
    }
    inputFile.close();

// Restarting Code


    // Compute the initial energy of the system
    double totalEnergy = CalculateTotalEnergy(xCoordinates, yCoordinates, zCoordinates, epsilon, sigma, cutoff, boxSize);
    std::cout << "Initial energy: " << totalEnergy << std::endl;

    // Monte Carlo simulation parameters
    int numSteps = 10000;
    double delta = 0.1; // Metropolis sampling width
    int rdfSamplingInterval = 10; 

    // Initialize random number generator
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> rand01(0.0, 1.0);
    
    std::ofstream energyFile("Total_Energy_V6.txt");
    std::ofstream coordinatesFile("Coordinates_V6.txt");
    std::ofstream rdfFile("rdf_values_averaged.txt");
    
    if (rdfFile.is_open()) {
        const double dr = 0.01; // Bin width
        const double maxRadius = a / 2;
        const int numBins = static_cast<int>(maxRadius / dr);
        std::vector<double> rdf(numBins, 0);

        // Perform the Monte Carlo simulation
        
        for (int step = 0; step < numSteps; ++step) {
            // Randomly select a particle
            int randomParticle = rand01(mt) * xCoordinates.size();
            double oldX = xCoordinates[randomParticle];
            double oldY = yCoordinates[randomParticle];
            double oldZ = zCoordinates[randomParticle];

            // Generate a random move
            double dx = (rand01(mt) - 0.5) * delta;
            double dy = (rand01(mt) - 0.5) * delta;
            double dz = (rand01(mt) - 0.5) * delta;

            // Apply periodic boundary conditions to the move
            double newX = oldX + dx;
            double newY = oldY + dy;
            double newZ = oldZ + dz;
            if (newX < 0) newX =newX + boxSize;
            if (newX > boxSize) newX =newX- boxSize;
            if (newY < 0) newY =newY+ boxSize;
            if (newY > boxSize) newY =newY- boxSize;
            if (newZ < 0) newZ =newZ + boxSize;
            if (newZ > boxSize) newZ = newZ - boxSize;


            // Calculate the energy change due to the move
            double newEnergy = totalEnergy - CalculateTotalEnergy(xCoordinates, yCoordinates, zCoordinates, epsilon, sigma, cutoff, boxSize);
            xCoordinates[randomParticle] = newX;
            yCoordinates[randomParticle] = newY;
            zCoordinates[randomParticle] = newZ;
            double energyChange = CalculateTotalEnergy(xCoordinates, yCoordinates, zCoordinates, epsilon, sigma, cutoff, boxSize) - newEnergy;

            // Accept or reject the move based on the Metropolis criterion
            if (energyChange < 0 || rand01(mt) < std::exp(-energyChange / temperature)) {
                totalEnergy += energyChange;
            } else {
                // Revert the move
                xCoordinates[randomParticle] = oldX;
                yCoordinates[randomParticle] = oldY;
                zCoordinates[randomParticle] = oldZ;
            }

            const double dr = 0.01; //bin width
            const double maxRadius = a/2; 
            const int numBins = static_cast<int>(maxRadius/ dr);
            const int NumCoords = xCoordinates.size();

            std::vector<double>rdf(numBins,0);
            std::vector<int> rdfCount(numBins, 0);

            if (step % rdfSamplingInterval == 0) {
                // Initialize RDF for this step
                std::vector<int> rdfStep(numBins, 0);

                // Calculating periodic boundary condtion distances
                std::vector<std::vector<double>> Distance(totalCoordinates, std::vector<double>(totalCoordinates, 0.0));

                for (int i = 0; i<xCoordinates.size(); ++i) {
                    for (int j = i+1; j<xCoordinates.size(); ++j) {
                        Distance[i][j] = PBC_Distance(xCoordinates[i], yCoordinates[i], zCoordinates[i], xCoordinates[j], yCoordinates[j], zCoordinates[j], a);
                    }
                }

                for (int i=0; i< NumCoords ; ++i){
                    for (int j=i+1; j< NumCoords ; ++j){
                        double distance = Distance[i][j];
                        if (distance <= maxRadius) {
                            int bin = static_cast<int>(distance / dr);
                            rdf[bin]++;
                        }
                    }
                }
                    for (int i=0;i<numBins;++i){
                        //std::cout<< i <<" "<<rdf[i] << " ,";
                        double r_out = (i+1) * dr;
                        double r_in = i*dr;
                        double OuterShellVolume = 4.0 /3 * M_PI * (r_out * r_out * r_out); // Outer Shell volume
                        double InnerShellVolume = 4.0 /3 * M_PI * (r_in * r_in * r_in); // Shell volume
                        //std::cout<< OuterShellVolume;
                        double shellvolume = OuterShellVolume - InnerShellVolume;
                        rdf[i] = rdf[i]/ (2700 *shellvolume);
                        //std::cout<< rdf[i]<<" ";
                    }
                // Add the RDF data for this step to the total
                //for (int i = 0; i < numBins; ++i) {
                //    rdf[i] += static_cast<double>(rdfStep[i]);
                //   rdfCount[i]++;
                //}
                if (step % rdfSamplingInterval == 0) {
            
                    for (int i = 0; i < numBins; ++i) {
                        double r = i * dr + dr / 2.0;
                        rdfFile << r << " " << rdf[i] << "\n";
                    }
                }
            } 
            rdfFile.close();

            // Record energy over moves
            if (step % 100 == 0) {
                std::cout << "Step " << step << ", Energy: " << totalEnergy << std::endl;
                energyFile << totalEnergy << std::endl;
                coordinatesFile << "2700 " << std::endl;
                coordinatesFile << "Step " << step << std::endl;
                for (int i = 0; i < xCoordinates.size(); ++i) {
                    coordinatesFile << "H  " << xCoordinates[i] << " " << yCoordinates[i] << " " << zCoordinates[i] << std::endl;
                }
            }
        }
    }
    energyFile.close();
    coordinatesFile.close();
    std::cout << "Final energy: " << totalEnergy << std::endl;

    return 0;
}}