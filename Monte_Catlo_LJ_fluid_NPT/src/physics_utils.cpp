#include <iostream>
#include <cmath>
#include <random>
#include <cstdlib>
#include <time.h>
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
    int i, j, k;
    const double Tr = 1.0;  // Reduced Temperature

    double Pr = 0.989;    // Reduced pressure
    static double Coord[N_Particles][3] = { 0 };
    srand(time(0));
    double N_density_new;
    double BoxVolume;
    
    // Read starting coordinates
    ifstream xyzFile("Atom_coordinates.txt");
    if (!xyzFile.is_open()) {
        cerr << "Error: Run setup script first to generate Atom_coordinates.txt" << endl;
        return 1;
    }
    
    int num_particles;
    xyzFile >> num_particles;
    string line;
    getline(xyzFile, line); // consume newline
    getline(xyzFile, line); // consume header
    for (i = 0; i < N_Particles; ++i)
    {
        char atom_type;
        xyzFile >> atom_type >> Coord[i][0] >> Coord[i][1] >> Coord[i][2];
    }
    xyzFile.close();

    // Computing Potential Energy for initial configuration
    double Initial_U = 0;
    for (i = 0; i < N_Particles - 1; i++){
        for (j = i + 1; j < N_Particles; j++){
            double dist = PBC_Distance(Coord[i][0], Coord[i][1], Coord[i][2],
                                       Coord[j][0], Coord[j][1], Coord[j][2], box_size);
            if (dist <= LJ_cutoff && dist != 0) {
                Initial_U += LJ_Potential(dist, 1.0, 1.0, LJ_cutoff);
            }
        }
    }
    
    // NPT Monte Carlo loop
    double RandomMovetype;
    double disp = 0.3;
    double vmax = 0.01;
    const int TotalMoves = 2000000;
    int Move;
    
    ofstream outFile2("Move_Energy_BoxSize_Density.txt");
    ofstream outFile3("Trajectory.txt");

    cout << "Starting NPT Simulation..." << endl;

    for (Move = 0; Move < TotalMoves; Move++) {
        static double Coord_1[3] = { 0 };
        double Dist_2, Dist_3;
        static double Uold, Unew;
        double P_acc;
        double Ediff = 0;
        double log_NewBoxVolume;
        double NewBoxVolume, Newbox_size;

        RandomMovetype = rand() % (N_Particles + 1) + 1;
        
        // Random translation move
        if (RandomMovetype < N_Particles) { 

            k = rand() % N_Particles;
            Coord_1[0] = Coord[k][0] + (((double)rand() / RAND_MAX) - 0.5) * disp;
            Coord_1[1] = Coord[k][1] + (((double)rand() / RAND_MAX) - 0.5) * disp;
            Coord_1[2] = Coord[k][2] + (((double)rand() / RAND_MAX) - 0.5) * disp;

            // Applying PBC
            Coord_1[0] = Coord_1[0] - box_size * floor((Coord_1[0]) / box_size);
            Coord_1[1] = Coord_1[1] - box_size * floor((Coord_1[1]) / box_size);
            Coord_1[2] = Coord_1[2] - box_size * floor((Coord_1[2]) / box_size);

            for (j = 0; j < N_Particles; j++) {
                Uold = 0;
                Unew = 0;
                if (j != k){
                    Dist_2 = PBC_Distance(Coord[j][0], Coord[j][1], Coord[j][2], 
                                          Coord[k][0], Coord[k][1], Coord[k][2], box_size);
                    if (Dist_2 <= LJ_cutoff) {
                        Uold = LJ_Potential(Dist_2, 1.0, 1.0, LJ_cutoff);
                    }
                    
                    Dist_3 = PBC_Distance(Coord[j][0], Coord[j][1], Coord[j][2], 
                                          Coord_1[0], Coord_1[1], Coord_1[2], box_size);
                    if (Dist_3 <= LJ_cutoff) {
                        Unew = LJ_Potential(Dist_3, 1.0, 1.0, LJ_cutoff);
                    }
                    Ediff += Unew - Uold;
                }
            }

            if (Ediff <= 0 || ((double)rand() / RAND_MAX) < exp(-Ediff / Tr)){
                Coord[k][0] = Coord_1[0];
                Coord[k][1] = Coord_1[1];
                Coord[k][2] = Coord_1[2];
                Initial_U += Ediff;
            } 
        } 
        // Volume change Monte Carlo move
        else {
            double scaling_factor;
            BoxVolume = pow(box_size, 3.0);
            log_NewBoxVolume = log(BoxVolume) + (((double)rand() / RAND_MAX) - 0.5) * vmax;
            NewBoxVolume = exp(log_NewBoxVolume);
            Newbox_size = pow(NewBoxVolume, 1.0 / 3.0);

            scaling_factor = Newbox_size / box_size;

            // BUG FIX: Correctly scale X, Y, and Z
            for (i = 0; i < N_Particles; i++) {
                Coord[i][0] *= scaling_factor;
                Coord[i][1] *= scaling_factor;
                Coord[i][2] *= scaling_factor;
            }

            double U_volchange = 0;
            
            for (i = 0; i < N_Particles - 1; i++){
                for (j = i + 1; j < N_Particles; j++){
                    double Dist_4 = PBC_Distance(Coord[i][0], Coord[i][1], Coord[i][2],
                                                 Coord[j][0], Coord[j][1], Coord[j][2], Newbox_size);
                    if (Dist_4 <= LJ_cutoff && Dist_4 != 0) {
                        U_volchange += LJ_Potential(Dist_4, 1.0, 1.0, LJ_cutoff);
                    }
                }
            }

            P_acc = -((U_volchange - Initial_U) + Pr * (NewBoxVolume - BoxVolume) - (N_Particles + 1) * log(scaling_factor));
            
            if ( exp(P_acc) < ((double)(rand() % 10 + 1) / 10.0) ){
                // BUG FIX: Rejecting volume change move requires dividing by scaling_factor to revert
                for (i = 0; i < N_Particles; i++){
                    Coord[i][0] /= scaling_factor;
                    Coord[i][1] /= scaling_factor;
                    Coord[i][2] /= scaling_factor;
                }
            } else {
                // Accepting volume change move
                Initial_U = U_volchange;
                box_size = Newbox_size;
                N_density_new = N_Particles / (pow(box_size, 3));
            }
        }

        // Recording thermodynamic data at every 100 steps
        if (Move % 100 == 0){
            outFile2 << Move << "   " << Initial_U << "   " <<  box_size << "   " << N_density_new << "\n";
        }
        
        // Recording trajectory data at every 1000 steps
        if (Move % 1000 == 0){
            outFile3 << N_Particles << "\n";
            outFile3 << "Generated XYZ coordinates\n";
            for (int m = 0; m < N_Particles; m++){
                outFile3 << "H" << " " << Coord[m][0] << "  " << Coord[m][1] << "  " << Coord[m][2] << endl;
            }
        }

        // Final RDF calculation at the end of equilibration
        if (Move == TotalMoves - 1){
            const double dr = 0.01; //bin width
            const double maxRadius = 8.0; 
            const int numBins = static_cast<int>(maxRadius / dr);

            vector<double> rdf(numBins, 0);
            
            for (i = 0; i < N_Particles - 1; i++ ){
                for(j = i + 1; j < N_Particles; j++){
                    double distance = PBC_Distance(Coord[i][0], Coord[i][1], Coord[i][2],
                                                   Coord[j][0], Coord[j][1], Coord[j][2], box_size);

                    if (distance <= maxRadius) {
                        int bin = static_cast<int>(distance / dr);
                        if (bin < numBins) rdf[bin] += 2;
                    }
                }
            }

            ofstream rdfFile("rdf_values.txt");
            if (rdfFile.is_open()) {
                for (int b = 0; b < numBins; ++b) {
                    double r_out = (b + 1) * dr;
                    double r_in = b * dr;
                    double shellvolume = (4.0 / 3.0) * M_PI * (pow(r_out, 3) - pow(r_in, 3));
                    rdf[b] = rdf[b] / (N_Particles * shellvolume * (N_Particles / pow(box_size, 3)));

                    double r = b * dr + dr / 2.0;
                    rdfFile << r << " " << rdf[b] << "\n";
                }
                rdfFile.close();
            }     
        }  
    } 
    outFile2.close();
    outFile3.close();
    cout << "Simulation Complete." << endl;

    return 0;
}