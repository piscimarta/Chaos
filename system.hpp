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
    //arma::vec compute_tot_spec_ang_mom();
    arma::vec compute_spec_ang_mom_2body();
    arma::vec compute_spec_ang_mom_3body(int j); //for spec. ang. mom with index j (m2->j=1); (m3->j=2)    
    double compute_eccentricity_2body();
    double compute_eccentricity_3body(int j); //computes ecc for m2(j=1) or m3 (j=2)
    double compute_semi_maj_ax_2body();
    double compute_semi_maj_ax_3body(int j); //computes sma for m2(j=1) or m3(j=2) 
    
    void evolveEuler(double h);
    void evolveEulerCromer(double h);
    void evolveLeapFrog(double h, int i, int N);
    void evolve_RK4(double h);
    
    void add_planet(double m, arma::vec r, arma::vec v);
    void initialize_kepler_orbit_2body(double e, double a, double m1, double m2);
    void initialize_kepler_orbit_3body(double e, double a1, double a2, double m1, double m2, double m3);
    double adaptive_time_step_diff_quot(double eta, double h);
};



#endif