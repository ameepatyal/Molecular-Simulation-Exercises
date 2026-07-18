#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include "../include/physics_utils.h"

using namespace std;

int main()
{
    double box_size = 18.0;
    double LJ_cutoff = 2.5; // LJ cutoff distance
    const int N_Particles = 4374; 
    int N_part = 0;  // Particle counter
    int x, y, z;

    static double Coord[N_Particles][3] = { 0 };
    const int num_bins = 100; // Number of bins for RDF
    double bin_width = LJ_cutoff / num_bins;
    double coords[N_Particles][3];

    // Using a simple cubic lattice to initialize coordinates
    // 4374 atoms on a x*y*z grid 16*17*17
    
    ofstream outFile("Atom_coordinates.txt");
    outFile << N_Particles << "\n";
    outFile << "Generated XYZ coordinates\n";

    for (x = 0; x < 17; x++){
        for (y = 0; y < 17; y++){
            for (z = 0; z < 16; z++){
                if (N_part < N_Particles){
                    Coord[N_part][0] = x * box_size / 17.0; 
                    Coord[N_part][1] = y * box_size / 17.0; 
                    Coord[N_part][2] = z * box_size / 16.0; 

                    outFile << "H" << " " << Coord[N_part][0] << "  " << Coord[N_part][1] << "  " << Coord[N_part][2] << endl;
                    ++N_part;
                }
            }
        }
    }
    outFile.close();
    cout << "Initial coordinates saved to 'Atom_coordinates.txt'" << endl;

    // Read coordinates back to verify and compute initial RDF
    ifstream xyzFile("Atom_coordinates.txt");
    int num_particles;
    xyzFile >> num_particles;
    string line;
    getline(xyzFile, line); // consume newline
    getline(xyzFile, line); // consume header

    for (int i = 0; i < N_Particles; ++i)
    {
        char atom_type;
        xyzFile >> atom_type >> coords[i][0] >> coords[i][1] >> coords[i][2];
    }
    xyzFile.close();

    vector<double> RDF(num_bins, 0.0);

    // Compute RDF
    for (int l = 0; l < N_Particles - 1; l++)
    {
        for (int t = l + 1; t < N_Particles; t++)
        {
            // Utilizing shared physics utility for distance
            double r_ij = PBC_Distance(coords[l][0], coords[l][1], coords[l][2], 
                                       coords[t][0], coords[t][1], coords[t][2], box_size);

            // Check if particle pair is within LJ cutoff (or max RDF radius)
            if (r_ij <= 8.0)
            {
                int bin_index = static_cast<int>(r_ij / bin_width);
                if (bin_index < num_bins) {
                    RDF[bin_index] += 2; // Count each pair twice for normalization
                }
            }
        }
    }

    // Normalize RDF
    for (int l = 0; l < num_bins; l++)
    {
        double r = l * bin_width + bin_width / 2.0;
        double shell_volume = 4.0 / 3.0 * M_PI * (pow(r + bin_width, 3) - pow(r, 3));
        RDF[l] /= (shell_volume * (N_Particles / pow(box_size, 3))); 
    }

    // Output RDF to file
    ofstream rdfFile("rdf_Initial_data.txt");
    for (int l = 0; l < num_bins; l++)
    {
        double r = l * bin_width + bin_width / 2.0;
        rdfFile << r << "\t" << RDF[l] << "\n";
    }
    rdfFile.close();
    cout << "Initial RDF saved to 'rdf_Initial_data.txt'" << endl;

    return 0;
}