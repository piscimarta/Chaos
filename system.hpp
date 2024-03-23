// include guard
#ifndef __system_hpp__   
#define __system_hpp__

// inbclude headers
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <armadillo>
#include "planet.hpp"

class System{
public:
    // variables 
    std::vector<Planet> planets;
    double G = 1;

    // constructors (?)

    // computing accelaration on the jth planet
    arma::vec compute_acceleration(int j);
    double compute_energy();
    void coord_transf(); 
    arma::vec compute_spec_ang_mom();
    double compute_eccentricity();
    double compute_semi_maj_ax();
    void evolveEuler(double h);
    void evolve_RK4(double h);
    void add_planet(double m, arma::vec r, arma::vec v);
    void initialize_kepler_orbit(double e, double a, double m1, double m2);
};



#endif