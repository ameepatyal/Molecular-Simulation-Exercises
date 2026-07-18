#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include "physics_utils.h"


int main() {

    double a = 15.0; // box size /edge length
    double nDensity = 0.8; // Particle number density
    int Gridintervals = 14; // Number of intervals
    double Grid_width = a / Gridintervals;
    int totalCoordinates = Gridintervals * Gridintervals * Gridintervals ;
    int numCoordinatesToDelete = totalCoordinates - (a*a*a) * nDensity +1; //44: We need to delete extra coordinates because 2744 particles lead to an increased density
    

    const double epsilon = 1.0; // Depth of the potential well
    const double sigma = 1.0;   // Finite distance at which potential is zero
    const double cutoff = 2.5 * sigma; // Cutoff distance

    //std::cout<<totalCoordinates << std::endl ;
    //std::cout<<numCoordinatesToDelete << std::endl ;

    std::vector<std::string> coordinates;

    std::vector<double> xCoordinates;
    std::vector<double> yCoordinates;
    std::vector<double> zCoordinates;

    for (int i = 0; i < Gridintervals; ++i) {
        for (int j = 0; j < Gridintervals; ++j) {
            for (int k = 0; k < Gridintervals; ++k) {
                double x = i * Grid_width;
                double y = j * Grid_width;
                double z = k * Grid_width;
                xCoordinates.push_back(x);
                yCoordinates.push_back(y);
                zCoordinates.push_back(z);
                coordinates.push_back("H " + std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z));
            }
        }
    }
   

    //Writing Initial cooridinates with 2744 atoms
    std::ofstream latticeFile("Initial_lattice_coordinates.txt");

    if (latticeFile.is_open()) {
        for (int i = 0; i < totalCoordinates; ++i) {
            latticeFile << coordinates[i] << std::endl;
        }

        latticeFile.close();
    }

    // Read coordinates from "Initial_lattice_coordinates.txt"
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(0, totalCoordinates - 1);
    
    std::ifstream inputFile("Initial_lattice_coordinates.txt");
    std::vector<std::string> inputCoordinates;
    std::string line;

    while (std::getline(inputFile, line)) {
        inputCoordinates.push_back(line);
    }

    inputFile.close();

    // Randomly delete coordinates to make number density 0.8
    
    std::vector<int> indicesToDelete;

    for (int i = 0; i < numCoordinatesToDelete; ++i) {
        int randomIndex = dis(gen);
        indicesToDelete.push_back(randomIndex);
    }

    for (int i = indicesToDelete.size() - 1; i >= 0; --i) {

        inputCoordinates.erase(inputCoordinates.begin() + indicesToDelete[i]);
        xCoordinates.erase(xCoordinates.begin()+ indicesToDelete[i]);
        yCoordinates.erase(yCoordinates.begin()+ indicesToDelete[i]);
        zCoordinates.erase(zCoordinates.begin()+ indicesToDelete[i]);
    }
    
    // Export to "Initial_NVT_coordinates.txt" without empty lines
    std::ofstream outputFile("Initial_NVT_coordinates.txt");

    if (outputFile.is_open()) {
        for (const std::string& coordinate : inputCoordinates) {
            outputFile << coordinate << std::endl;
        }

        outputFile.close();
    }

    // Calculating periodic boundary condtion distances
    std::vector<std::vector<double>> Distance(totalCoordinates, std::vector<double>(totalCoordinates, 0.0));

    for (int i = 0; i<xCoordinates.size(); ++i) {
        for (int j = i+1; j<xCoordinates.size(); ++j) {
            Distance[i][j] = PBC_Distance(xCoordinates[i], yCoordinates[i], zCoordinates[i], xCoordinates[j], yCoordinates[j], zCoordinates[j], a);
        }
    }
    
    //Radial Distribution Function

    const double dr = 0.01; //bin width
    const double maxRadius = a/2; 
    const int numBins = static_cast<int>(maxRadius/ dr);
    const int NumCoords = xCoordinates.size();

    std::vector<double>rdf(numBins,0);

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
    
    std::ofstream rdfFile("rdf_values.txt");
    if (rdfFile.is_open()) {
        for (int i = 0; i < numBins; ++i) {
            double r = i * dr + dr / 2.0;
            rdfFile << r << " " << rdf[i] << "\n";
        }
        rdfFile.close();
    }

    // Compute total energy
    double totalEnergy = 0.0;
    for (int i = 0; i < xCoordinates.size(); ++i) {
        for (int j = i + 1; j < xCoordinates.size(); ++j) {
            double distance = Distance[i][j];
            if (distance <= maxRadius) {
                totalEnergy += LJPotential(distance, epsilon,  sigma,  cutoff);
            }
        }
    }

    std::cout << "Total Energy: " << totalEnergy << std::endl;
    

    return 0;
}