// include header
#include "planet.hpp"
#include "system.hpp"
#include "header.hpp"

// useful line for debugging:
//std::cout << __FILE__ << __LINE__<< std::endl;

int main(){

// 
double h = 2*M_PI/10000;
int iter = 100000;

// planet parameters 
double m1 = 1;
double m2 = 1e-3;
double e = 0.5;
double a = 1;

// printing parameters
double width = 7;
double prec = 7;

// create a system Sosy
System Sosy;

Sosy.initialize_kepler_orbit(e, a, m1 ,m2);
Sosy.coord_transf();

// open file .txt to save data
std::ofstream ofile;
ofile.open("test3.txt");

// save to file data: t p1 x_1 y_1 z_1 vx_1 vy_1 vz_1 p2 x_2 y_2 z_2 vx_2 vy_2 vz_2
for(int i = 0; i < iter ; i++){
    ofile  << scientific_format(h*i, width, prec);

    for(int i_planet = 0 ; i_planet < Sosy.planets.size(); i_planet++){
    ofile   << " " << i_planet+1
                    << " " << scientific_format(Sosy.planets.at(i_planet).r(0), width, prec)
                    << " " << scientific_format(Sosy.planets.at(i_planet).r(1), width, prec)
                    << " " << scientific_format(Sosy.planets.at(i_planet).r(2), width, prec)
                    << " " << scientific_format(Sosy.planets.at(i_planet).v(0), width, prec)
                    << " " << scientific_format(Sosy.planets.at(i_planet).v(1), width, prec)
                    << " " << scientific_format(Sosy.planets.at(i_planet).v(2), width, prec);

    }
    ofile << std::endl;

   // Evolve the system with Euler 
   Sosy.evolveEuler(h);
    
}
    // close file
    ofile.close();
    return 0;
}