#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <chrono>

using namespace std;
using dur = std::chrono::high_resolution_clock::duration;

// ------------------------------------------------------------
// Global parameters
// ------------------------------------------------------------
const double gamma_val = 1.4;   // Ratio of specific heats
const double CFL = 0.5;         // CFL number

// ------------------------------------------------------------
// Compute pressure from the conservative variables
// ------------------------------------------------------------
double pressure(double rho, double rhou, double rhov, double E) {
    double u = rhou / rho;
    double v = rhov / rho;
    double kinetic = 0.5 * rho * (u * u + v * v);
    return (gamma_val - 1.0) * (E - kinetic);
}

// ------------------------------------------------------------
// Compute flux in the x-direction
// ------------------------------------------------------------
void fluxX(double rho, double rhou, double rhov, double E, 
           double& frho, double& frhou, double& frhov, double& fE) {
    double u = rhou / rho;
    double p = pressure(rho, rhou, rhov, E);
    frho = rhou;
    frhou = rhou * u + p;
    frhov = rhov * u;
    fE = (E + p) * u;
}

// ------------------------------------------------------------
// Compute flux in the y-direction
// ------------------------------------------------------------
void fluxY(double rho, double rhou, double rhov, double E,
           double& frho, double& frhou, double& frhov, double& fE) {
    double v = rhov / rho;
    double p = pressure(rho, rhou, rhov, E);
    frho = rhov;
    frhou = rhou * v;
    frhov = rhov * v + p;
    fE = (E + p) * v;
}

// ------------------------------------------------------------
// Main simulation routine
// ------------------------------------------------------------
int main(int argc, char* argv[]){
    // ----- Grid and domain parameters -----

	int mult = 1;

	if(argc == 2)
	{
	try
	{
		std::vector<std::string> args(argv,argv+argc);
		mult = std::stoi(args.at(1));
	}
	catch(...){}
	}

	std::cout << "Running at x" << mult << " size" << std::endl;

    const int Nx = 200*mult;         // Number of cells in x (excluding ghost cells)
    const int Ny = 100*mult;         // Number of cells in y
    const double Lx = 2.0;      // Domain length in x
    const double Ly = 1.0;      // Domain length in y
    const double dx = Lx / Nx;
    const double dy = Ly / Ny;

	//timing variables
	dur main_time(0);

	dur cupdate_time(0);
	int cupdate_count = 0;
	dur kinetic_time(0);
	int kinetic_count = 0;
	dur copy_time(0);
	int copy_count = 0;
	dur bound_time(0);
	int bound_count = 0;

    // Create flat arrays (with ghost cells)
    const int total_size = (Nx + 2) * (Ny + 2);
    
    vector<double> rho(total_size);
    vector<double> rhou(total_size);
    vector<double> rhov(total_size);
    vector<double> E(total_size);
    
    vector<double> rho_new(total_size);
    vector<double> rhou_new(total_size);
    vector<double> rhov_new(total_size);
    vector<double> E_new(total_size);
    
    // A mask to mark solid cells (inside the cylinder)
    vector<bool> solid(total_size, false);

    // ----- Obstacle (cylinder) parameters -----
    const double cx = 0.5;      // Cylinder center x
    const double cy = 0.5;      // Cylinder center y
    const double radius = 0.1;  // Cylinder radius

    // ----- Free-stream initial conditions (inflow) -----
    const double rho0 = 1.0;
    const double u0 = 1.0;
    const double v0 = 0.0;
    const double p0 = 1.0;
    const double E0 = p0/(gamma_val - 1.0) + 0.5*rho0*(u0*u0 + v0*v0);

	//this rund once and would take no time, so I'll skip it
    // ----- Initialize grid and obstacle mask -----
    for (int i = 0; i < Nx+2; i++){
        for (int j = 0; j < Ny+2; j++){
            // Compute cell center coordinates
            double x = (i - 0.5) * dx;
            double y = (j - 0.5) * dy;
            // Mark cell as solid if inside the cylinder
            if ((x - cx)*(x - cx) + (y - cy)*(y - cy) <= radius * radius) {
                solid[i*(Ny+2)+j] = true;
                // For a wall, we set zero velocity
                rho[i*(Ny+2)+j] = rho0;
                rhou[i*(Ny+2)+j] = 0.0;
                rhov[i*(Ny+2)+j] = 0.0;
                E[i*(Ny+2)+j] = p0/(gamma_val - 1.0);
            } else {
                solid[i*(Ny+2)+j] = false;
                rho[i*(Ny+2)+j] = rho0;
                rhou[i*(Ny+2)+j] = rho0 * u0;
                rhov[i*(Ny+2)+j] = rho0 * v0;
                E[i*(Ny+2)+j] = E0;
            }
        }
    }

    // ----- Determine time step from CFL condition -----
    double c0 = sqrt(gamma_val * p0 / rho0);
    double dt = CFL * min(dx, dy) / (fabs(u0) + c0)/2.0;

    // ----- Time stepping parameters -----
    const int nSteps = 2000;

    // ----- Main time-stepping loop -----
	auto main_start = std::chrono::high_resolution_clock::now();
    for (int n = 0; n < nSteps; n++){

		auto bound_start = std::chrono::high_resolution_clock::now();
        // --- Apply boundary conditions on ghost cells ---
        // Left boundary (inflow): fixed free-stream state
        for (int j = 0; j < Ny+2; j++){
            rho[0*(Ny+2)+j] = rho0;
            rhou[0*(Ny+2)+j] = rho0*u0;
            rhov[0*(Ny+2)+j] = rho0*v0;
            E[0*(Ny+2)+j] = E0;
			bound_count++;
        }
        // Right boundary (outflow): copy from the interior
        for (int j = 0; j < Ny+2; j++){
            rho[(Nx+1)*(Ny+2)+j] = rho[Nx*(Ny+2)+j];
            rhou[(Nx+1)*(Ny+2)+j] = rhou[Nx*(Ny+2)+j];
            rhov[(Nx+1)*(Ny+2)+j] = rhov[Nx*(Ny+2)+j];
            E[(Nx+1)*(Ny+2)+j] = E[Nx*(Ny+2)+j];
			bound_count++;
        }
        // Bottom boundary: reflective
        for (int i = 0; i < Nx+2; i++){
            rho[i*(Ny+2)+0] = rho[i*(Ny+2)+1];
            rhou[i*(Ny+2)+0] = rhou[i*(Ny+2)+1];
            rhov[i*(Ny+2)+0] = -rhov[i*(Ny+2)+1];
            E[i*(Ny+2)+0] = E[i*(Ny+2)+1];
			bound_count++;
        }
        // Top boundary: reflective
        for (int i = 0; i < Nx+2; i++){
            rho[i*(Ny+2)+(Ny+1)] = rho[i*(Ny+2)+Ny];
            rhou[i*(Ny+2)+(Ny+1)] = rhou[i*(Ny+2)+Ny];
            rhov[i*(Ny+2)+(Ny+1)] = -rhov[i*(Ny+2)+Ny];
            E[i*(Ny+2)+(Ny+1)] = E[i*(Ny+2)+Ny];
			bound_count++;
        }
		auto bound_end = std::chrono::high_resolution_clock::now();
		bound_time += bound_end-bound_start;

        // --- Update interior cells using a Lax-Friedrichs scheme ---
		auto cupdate_start = std::chrono::high_resolution_clock::now();
        
		//notice that each iteration only modifies their own  cell
		#pragma omp parallel for shared(solid,rho,rhou,rhov,E,rho_new,rhou_new,rhov_new,E_new) reduction(+:cupdate_count) collapse(2)
		for (int i = 1; i <= Nx; i++){
            for (int j = 1; j <= Ny; j++){
                // If the cell is inside the solid obstacle, do not update it
                if (solid[i*(Ny+2)+j]) {
                    rho_new[i*(Ny+2)+j] = rho[i*(Ny+2)+j];
                    rhou_new[i*(Ny+2)+j] = rhou[i*(Ny+2)+j];
                    rhov_new[i*(Ny+2)+j] = rhov[i*(Ny+2)+j];
                    E_new[i*(Ny+2)+j] = E[i*(Ny+2)+j];
                    continue;
                }

                // Compute a Lax averaging of the four neighboring cells
                rho_new[i*(Ny+2)+j] = 0.25 * (rho[(i+1)*(Ny+2)+j] + rho[(i-1)*(Ny+2)+j] + 
                                             rho[i*(Ny+2)+(j+1)] + rho[i*(Ny+2)+(j-1)]);
                rhou_new[i*(Ny+2)+j] = 0.25 * (rhou[(i+1)*(Ny+2)+j] + rhou[(i-1)*(Ny+2)+j] + 
                                              rhou[i*(Ny+2)+(j+1)] + rhou[i*(Ny+2)+(j-1)]);
                rhov_new[i*(Ny+2)+j] = 0.25 * (rhov[(i+1)*(Ny+2)+j] + rhov[(i-1)*(Ny+2)+j] + 
                                              rhov[i*(Ny+2)+(j+1)] + rhov[i*(Ny+2)+(j-1)]);
                E_new[i*(Ny+2)+j] = 0.25 * (E[(i+1)*(Ny+2)+j] + E[(i-1)*(Ny+2)+j] + 
                                           E[i*(Ny+2)+(j+1)] + E[i*(Ny+2)+(j-1)]);

                // Compute fluxes
                double fx_rho1, fx_rhou1, fx_rhov1, fx_E1;
                double fx_rho2, fx_rhou2, fx_rhov2, fx_E2;
                double fy_rho1, fy_rhou1, fy_rhov1, fy_E1;
                double fy_rho2, fy_rhou2, fy_rhov2, fy_E2;

                fluxX(rho[(i+1)*(Ny+2)+j], rhou[(i+1)*(Ny+2)+j], rhov[(i+1)*(Ny+2)+j], E[(i+1)*(Ny+2)+j],
                      fx_rho1, fx_rhou1, fx_rhov1, fx_E1);
                fluxX(rho[(i-1)*(Ny+2)+j], rhou[(i-1)*(Ny+2)+j], rhov[(i-1)*(Ny+2)+j], E[(i-1)*(Ny+2)+j],
                      fx_rho2, fx_rhou2, fx_rhov2, fx_E2);
                fluxY(rho[i*(Ny+2)+(j+1)], rhou[i*(Ny+2)+(j+1)], rhov[i*(Ny+2)+(j+1)], E[i*(Ny+2)+(j+1)],
                      fy_rho1, fy_rhou1, fy_rhov1, fy_E1);
                fluxY(rho[i*(Ny+2)+(j-1)], rhou[i*(Ny+2)+(j-1)], rhov[i*(Ny+2)+(j-1)], E[i*(Ny+2)+(j-1)],
                      fy_rho2, fy_rhou2, fy_rhov2, fy_E2);

                // Apply flux differences
                double dtdx = dt / (2 * dx);
                double dtdy = dt / (2 * dy);
                
                rho_new[i*(Ny+2)+j] -= dtdx * (fx_rho1 - fx_rho2) + dtdy * (fy_rho1 - fy_rho2);
                rhou_new[i*(Ny+2)+j] -= dtdx * (fx_rhou1 - fx_rhou2) + dtdy * (fy_rhou1 - fy_rhou2);
                rhov_new[i*(Ny+2)+j] -= dtdx * (fx_rhov1 - fx_rhov2) + dtdy * (fy_rhov1 - fy_rhov2);
                E_new[i*(Ny+2)+j] -= dtdx * (fx_E1 - fx_E2) + dtdy * (fy_E1 - fy_E2);

				cupdate_count++;
            }
        }
		auto cupdate_end = std::chrono::high_resolution_clock::now();
		cupdate_time += cupdate_end-cupdate_start;

        // Copy updated values back
		auto copy_start = std::chrono::high_resolution_clock::now();
		#pragma omp parallel for shared(rho,rhou,rhov,E,rho_new,rhou_new,rhov_new,E_new) reduction(+:copy_count) collapse(2)
        for (int i = 1; i <= Nx; i++){
            for (int j = 1; j <= Ny; j++){
                rho[i*(Ny+2)+j] = rho_new[i*(Ny+2)+j];
                rhou[i*(Ny+2)+j] = rhou_new[i*(Ny+2)+j];
                rhov[i*(Ny+2)+j] = rhov_new[i*(Ny+2)+j];
                E[i*(Ny+2)+j] = E_new[i*(Ny+2)+j];
				copy_count++;
            }
        }
		auto copy_end = std::chrono::high_resolution_clock::now();
		copy_time += copy_end-copy_start;

        // Calculate total kinetic energy
        double total_kinetic = 0.0;
		auto kinetic_start = std::chrono::high_resolution_clock::now();
		#pragma omp parallel for shared(rho,rhou,rhov) reduction(+:kinetic_count,total_kinetic) collapse(2)
        for (int i = 1; i <= Nx; i++) {
            for (int j = 1; j <= Ny; j++) {

                double u = rhou[i*(Ny+2)+j] / rho[i*(Ny+2)+j];
                double v = rhov[i*(Ny+2)+j] / rho[i*(Ny+2)+j];
                total_kinetic += 0.5 * rho[i*(Ny+2)+j] * (u * u + v * v);
				kinetic_count++;
            }
        }
		auto kinetic_end = std::chrono::high_resolution_clock::now();
		kinetic_time += kinetic_end-kinetic_start;

        if (n % 50 == 0) {
            cout << "Step " << n << " completed, total kinetic energy: " << total_kinetic << endl;
        }


    }
	auto main_end = std::chrono::high_resolution_clock::now();
	main_time += main_end-main_start;

	std::cout << "Execution time of main time step: " << std::chrono::duration_cast<std::chrono::milliseconds>(main_time).count() << " ms"<< std::endl;
	std::cout << "Execution time of boundary conditions: " << std::chrono::duration_cast<std::chrono::milliseconds>(bound_time).count() << " ms"<< std::endl;
	std::cout << "Execution time of cell updates: " << std::chrono::duration_cast<std::chrono::milliseconds>(cupdate_time).count() << " ms" <<  std::endl;
	std::cout << "Execution time of update copies: " << std::chrono::duration_cast<std::chrono::milliseconds>(copy_time).count() << " ms" << std::endl;
	std::cout << "Execution time of kinetic energy calculation: " << std::chrono::duration_cast<std::chrono::milliseconds>(kinetic_time).count() << " ms" << std::endl;

    return 0;
}

