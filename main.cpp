#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>

// Simulation parameters
const int Nx = 400;    // resolution x-dir
const int Ny = 100;    // resolution y-dir
const int NL = 9;
const int Nt = 4000;   // number of timesteps
const double rho0 = 100.0;   // average density
const double tau = 0.6;      // collision timescale

// Lattice speeds / weights
int idxs[NL] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
int cxs[NL] = {0, 0, 1, 1, 1, 0, -1, -1, -1};
int cys[NL] = {0, 1, 1, 0, -1, -1, -1, 0, 1};
double weights[NL] = {4.0/9, 1.0/9, 1.0/36, 1.0/9, 1.0/36, 1.0/9, 1.0/36, 1.0/9, 1.0/36};

// Helper functions
int index(int i, int j, int k) {
    return i * Ny * NL + j * NL + k;
}

int index_2d(int i, int j) {
    return i * Ny + j;
}

// Main simulation function
void simulate() {
    std::vector<double> F(Nx * Ny * NL, 1.0);  // Initialize with ones
    std::vector<double> rho(Nx * Ny, 0.0);
    std::vector<double> Feq(Nx * Ny * NL, 0.0); // Initialize Feq

    // Initial Conditions - flow to the right with some perturbations
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            F[index(i, j, 3)] += 2.0 * (1.0 + 0.2 * std::cos(2 * M_PI * i / Nx * 4));
            for (int k = 0; k < NL; ++k) {
                F[index(i, j, k)] *= rho0 / rho[index(i, j, k)];
            }
        }
    }

    // Cylinder boundary
    std::vector<bool> cylinder(Nx * Ny, false);
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            cylinder[index_2d(i, j)] = (std::pow(i - Nx / 4, 2) + std::pow(j - Ny / 2, 2) < std::pow(Ny / 4, 2));
        }
    }


    // Simulation Main Loop
    for (int it = 0; it < Nt; ++it) {
        // Drift
        for (int k = 0; k < NL; ++k) {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    F[index(i, j, k)] = F[index((i - cxs[k] + Nx) % Nx, (j - cys[k] + Ny) % Ny, k)];
                }
            }
        }

        // Set reflective boundaries
        std::vector<double> bndryF;
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                if (cylinder[index_2d(i, j)]) {
                    {
                        for (int k = 0; k < NL; ++k) {
                            bndryF.push_back(F[index(i, j, (k + 3) % NL)]);
                        }
                    }
                }
            }
        }

        // Calculate fluid variables
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                rho[index_2d(i, j)] = 0.0;
                for (int k = 0; k < NL; ++k) {
                    rho[index_2d(i, j)] += F[index(i, j, k)];
                }
            }
        }

        // Apply Collision
        for (int k = 0; k < NL; ++k) {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    Feq[index(i, j, k)] = rho[index_2d(i, j)] * weights[k] *
                                          (1 + 3 * (cxs[k] + cys[k]) + 9 * std::pow((cxs[k] + cys[k]), 2) / 2 -
                                           3 * (std::pow(cxs[k], 2) + std::pow(cys[k], 2)) / 2);
                }
            }
        }

        for (int k = 0; k < NL; ++k) {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    F[index(i, j, k)] += -(1.0 / tau) * (F[index(i, j, k)] - Feq[index(i, j, k)]);
                }
            }
        }

        // Apply boundary
        int bndryIndex = 0;
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                if (cylinder[index_2d(i, j)]) {
                    {
                        for (int k = 0; k < NL; ++k) {
                            F[index(i, j, k)] = bndryF[bndryIndex++];
                        }
                    }
                }
            }
        }
    }
}

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    simulate();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    std::cout << "Time elapsed for simulation: " << duration.count() << " seconds" << std::endl;

    return 0;
}