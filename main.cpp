
// include header
#include "planet.hpp"
#include "system.hpp"
#include "header.hpp"

int main(){

double h = 2*M_PI/10000;
System Sosy;
double m1 = 1;
double m2 = 1e-3;
double e = 0.5;
double a = 1;
int iter = 100000;

//printing parameters
double width = 7;
double prec = 7;


//std::cout << Sosy.planets.size() << std::endl;
Sosy.initialize_kepler_orbit(e, a, m1 ,m2);
Sosy.coord_transf();


//std::cout << Sosy.planets.size()<< std::endl;
std::ofstream ofile;
ofile.open("test3.txt");


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

    
   Sosy.evolveEuler(h);
    
}
    // close file
    ofile.close();
    return 0;
}