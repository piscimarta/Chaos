#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <cassert>
#include <armadillo>
#define _USE_MATH_DEFINES

// include headers
#include "system.hpp"
#include "planet.hpp"


 arma::vec System::compute_acceleration(int i){
    // 
    arma::vec a = arma::vec(3).fill(0.);
    arma::vec r_i = planets.at(i).r;

    for (int j = 0 ; j < planets.size() ; j++){
        // def rij
        arma::vec rij = planets.at(j).r - planets.at(i).r;
        // norm of rij
        double r = norm(rij);
        
        if (i==j){
        a += arma::vec(3).fill(0.);
    }
        else{
        a += G*planets.at(j).m*(rij/(r*r*r));
    }
    }

    return a; 
 }

 double System::compute_energy(){
    double E_kin = 0;
    double E_pot = 0;

    for (int i = 0 ; i < planets.size() ; i++){
        double v_norm = norm(planets.at(i).v);

        E_kin += 0.5*planets.at(i).m*v_norm*v_norm;

        for (int j = 0 ; j < planets.size() ; j++){
            // def rij
            arma::vec rij = planets.at(j).r - planets.at(i).r;
            // norm of rij
            double r = norm(rij);
        if (i==j){
             E_pot += 0;
    }
        else{
             E_pot += G*planets.at(j).m*planets.at(i).m/(r);
    }
    }
    }
    return E_kin + E_pot;
 }

// shifts the coordinate system into the CM
void System::coord_transf(){
    arma::vec r_cm = arma::vec(3).fill(0.);
    arma::vec v_cm = arma::vec(3).fill(0.);
    arma::vec r_new = arma::vec(3);
    arma::vec v_new = arma::vec(3);
    double M = 0; 

    for (int j = 0 ; j < planets.size() ; j++){
        M += planets.at(j).m;
    } 
    for (int j = 0 ; j < planets.size() ; j++){
        r_cm +=  (planets.at(j).m*planets.at(j).r)/M;
        v_cm +=  (planets.at(j).m*planets.at(j).v)/M;
    } 
    for (int j = 0 ; j < planets.size() ; j++){
        r_new = planets.at(j).r - r_cm;
        v_new = planets.at(j).v - v_cm;
        planets.at(j).r = r_new;
        planets.at(j).v = v_new;
   }
}

 arma::vec System::compute_spec_ang_mom(){
    arma::vec L = arma::vec(3).fill(0.);

    for (int i = 0 ; i < planets.size() ; i++){
        L += arma::cross(planets.at(i).r, planets.at(i).v);
    }
    return L;
 }
// Evolve the system by one time step, h, using Euler
void System::evolveEuler(double h){
    std::vector<Planet> evolved_planets; 
    arma::vec r_new = arma::vec(3);
    arma::vec v_new = arma::vec(3);
    Planet p; 

    for (int j = 0; j < planets.size() ; j++){
        // dr/dt = v ==> r_j+i = r_j + h*v_ij
        r_new = planets.at(j).r + planets.at(j).v * h;
        // dv/dt = a  ==>  v_j+1 = v_j + h*a_j
        v_new = planets.at(j).v + h * compute_acceleration(j);
        // save changes in new state of the particle
        p.r = r_new;
        p.v = v_new;
        p.m = planets.at(j).m;
         // append the particle at t+dt at the end of the new vector
        evolved_planets.push_back(p);
    }
    planets = evolved_planets;
}
// add a planet to the system specifying r, v, m 
void System::add_planet(double m, arma::vec r, arma::vec v){
     Planet p; 
     p.r = r;
     p.v = v;
     p.m = m;
     planets.push_back(p);
}
// set initial conditions of two planets in  kepler orbit
void System::initialize_kepler_orbit(double e, double a, double m1, double m2){
        arma::vec r = arma::vec(3);
        arma::vec v = arma::vec(3);
        add_planet(m1, r, v);
        add_planet(m2, r, v);
        planets.at(0).r = arma::vec(3).fill(0.);
        planets.at(0).v = arma::vec(3).fill(0.);
        planets.at(1).r(0) = a*(1+e);
        planets.at(1).r(1) = 0;
        planets.at(1).r(2) = 0;
        planets.at(1).v(0) = 0;
        planets.at(1).v(1) = sqrt((G*(m1 + m2)/a)*((1-e)/(1+e)));
        planets.at(1).v(2) = 0;
    }