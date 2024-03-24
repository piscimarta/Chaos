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

        for (int j = i+1 ; j < planets.size() ; j++){
            // def rij
            arma::vec rij = planets.at(j).r - planets.at(i).r;
            // norm of rij
            double r = norm(rij);
            if (i==j){
                E_pot += 0;
            }
            else{
                E_pot -= G*planets.at(j).m*planets.at(i).m/(r);
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

 double System::compute_semi_maj_ax(){
    double a = 0;
    double M = 0;
    double j = 0;
    double e = 0;
    M = planets.at(0).m - planets.at(1).m;
    j = norm(compute_spec_ang_mom());
    e = compute_eccentricity();
    a = (j*j/(G*M))/(1-e*e);
     return a;
}

 double System::compute_eccentricity(){
    arma::vec e = arma::vec(3).fill(0.);
    // def rij
    arma::vec rij = planets.at(0).r - planets.at(1).r;
    // def vij
    arma::vec vij = planets.at(0).v - planets.at(1).v;
    // norm of rij
    double r = norm(rij);
    e = (arma::cross(vij, compute_spec_ang_mom()))/(G*(planets.at(0).m+planets.at(1).m)) - (rij)/r;
    return norm(e);
 }



void System::evolveEuler(double h){
    std::vector<Planet> evolved_planets;
    Planet p; 
    arma::vec r_new = arma::vec(3);
    arma::vec v_new = arma::vec(3);
        arma::vec a = arma::vec(3);
        
    //save accel. so we do not need to calculate it twice 

    for (int j = 0; j < planets.size() ; j++){
        // dv/dt = a  ==>  v_j+1 = v_j + h*a_j
        a= compute_acceleration(j);
        v_new = planets.at(j).v + h * a;
        // dr/dt = v ==> r_j+1 = r_j + h*v_ij
        r_new = planets.at(j).r + planets.at(j).v * h;
        p.a_prev = a;
        //save the previous accel. before updating the positions
        //we need prev. acc to access it for the adaptable time_step, 
        //implement this for all the other evolve.. aswell!!!!!!!
        p.r = r_new;
        p.v = v_new;
        p.m = planets.at(j).m;
        evolved_planets.push_back(p);
    }
    planets = evolved_planets;
    //save the system in evolved planets bc we need the previous positions for a
    //then later update all planets at the same time

}




void System::evolveEulerCromer(double h){
    std::vector<Planet> evolved_planets;
    Planet p; 
    arma::vec r_new = arma::vec(3);
    arma::vec v_new = arma::vec(3);
    arma::vec a = arma::vec(3);
    //Planet p; 

    for (int j = 0; j < planets.size() ; j++){
        // dv/dt = a  ==>  v_j+1 = v_j + h*a_j
        a= compute_acceleration(j);
        v_new = planets.at(j).v + h * a;
        // dr/dt = v ==> r_j+1 = r_j + h*v_ij
        p.v = v_new;
        r_new = planets.at(j).r + v_new * h;
        p.a_prev = a;
        //save the previous accel. before updating the positions
        //we need prev. acc to access it for the adaptable time_step, 

        p.r = r_new;
        p.m = planets.at(j).m;
        evolved_planets.push_back(p);
        //only difference to Euler: v is updated before we calculate r
    }
    planets= evolved_planets;
}


void System::evolveLeapFrog(double h, int i, int N){
    //now we need the step i and the #total steps N bc the integration step depends on them
    std::vector<Planet> evolved_planets;
    Planet p;
    arma::vec r_new = arma::vec(3);
    arma::vec v_new =arma::vec(3);
    arma::vec a = arma::vec(3);

    //here we need separate loops for r & v bc we update them separately

    if (i==0){
        //loop for updating r & a_prev
        for(int j = 0; j<planets.size(); j++){
            a = compute_acceleration(j);
            p.a_prev=a; 
            //save a_0 as prev with actual values
            r_new = planets.at(j).r + planets.at(j).v * h/2;
            //r_1/2 = r_0 + v_0 *h/2        
            //r_1/2 = r[1] in array
            //update pos. first
            p.r = r_new;
            p.v = planets.at(j).v;
            p.m = planets.at(j).m;
            evolved_planets.push_back(p);
            //need to save the updated planets in intermediate vector bc otherwise it would change a_prev
        }
        planets =evolved_planets;

        //loop for updating v
        for(int j=0; j<planets.size(); j++){
            // calc a_1/2 = a (t_1/2, r_1/2)
            a= compute_acceleration(j);
            v_new = planets.at(j).v +a*h;
            //v_1 = v_0 + a_(1/2) *h
            //now calculate velocities with updated coordinates
            //here we also do not need intermediate steps
            planets.at(j).v =v_new;
        }
    }
        
    else if (i>0 && i<N-1){
        //loop for r and a_prev
        for(int j =0; j<planets.size(); j++){
            a= compute_acceleration(j);
            p.a_prev = a;
            //save the previous acc. before updating r
            r_new = planets.at(j).r + planets.at(j).v *h;
            //r_(i-1/2) = r_(i+1-1/2) + v_(i+1)*h
            //r_(i-1/2) = r[i] in array
            p.r = r_new;
            p.v = planets.at(j).v;
            p.m = planets.at(j).m;
            evolved_planets.push_back(p);
            //need to save the updated planets in intermediate vector bc otherwise it would change a_prev
        }
        planets = evolved_planets;
        //again: update all pos. first, then accel.->velocity
        //loop for v
        for(int j =0; j<planets.size(); j++){
            a = compute_acceleration(j); //a_(i+1/2)
            v_new = planets.at(j).v + a *h;
            //v_i = v_(i-1)         + a_(i-1/2) *h
            //v_i =v[i] in array
            planets.at(j).v = v_new;    
        }   
    }
        
    else{  // i=N-1, last step, only pos update
    //this condition is difficult to do with an adaptable timestep
    //to calculate the timestep we need to calculate the acc
        for(int j=0; j<planets.size(); j++){
            r_new = planets.at(j).r + planets.at(j).v *h/2;
            //  //r_(n+1) = r_(n+1/2) + v_(n+1) * h/2
            planets.at(j).r =r_new;
        }
    }
    //planets = evolved_planets;
}


        //Iter = #Iterations =N
        //i goes from 0 to N-1

        //r              a           v  
        //r_0            -           v_0
        //r_1/2          a_1/2       -
        //-             -            v_1
        //r_3/2         a_3/2        -
        //-             -            v_2
        //...
        // -            -            v_(N-2)
        //r_(N-2 +1/2)  a(N-2 +1/2)  -
        //r_(N-1)       -            v_(N-1)

        
        //r_(N-1) = r_(N-2 +1/2) + v_(N-2) *h/2
        //v_(N-2) = v_(N-3) + a_(N-3 +1/2) *h
        //v_(N-1) = v_(N-2) + a_(N-2 +1/2) *h



        //r: r_0, r_1/2, r_(1+1/2), ... r_(N-2 +1/2), r_(N-1)     #N+1 locations 
        //a:      a_1/2, a_(1+1/2), ... a_(N-2 +1/2)              #N-1 acc.
        //v: v_0,   v_1,       v_2, ... v_(N-1)                   #N   velocities
        //#:  0      1       2            N-1            N

        //accel. don't need to be saved 


void System::evolve_RK4(double h){
        std::vector<arma::vec> K1_r_all(planets.size());
        std::vector<arma::vec> K1_v_all(planets.size());
        std::vector<arma::vec> K2_r_all(planets.size());
        std::vector<arma::vec> K2_v_all(planets.size());
        std::vector<arma::vec> K3_r_all(planets.size());
        std::vector<arma::vec> K3_v_all(planets.size());
        std::vector<arma::vec> K4_r_all(planets.size());
        std::vector<arma::vec> K4_v_all(planets.size());

        std::vector<Planet> initial_state = planets;
        std::vector<Planet> final_state = planets;
        //k1
        for (int i = 0; i < planets.size() ; i++){
            K1_r_all.at(i) = planets.at(i).v;
            K1_v_all.at(i) = compute_acceleration(i);
        }

        //k2
        for(int i = 0; i < planets.size(); i++){
            planets.at(i).r += 0.5*h*K1_r_all.at(i); 
            planets.at(i).v += 0.5*h*K1_v_all.at(i);
        }
        //evolveEuler(0.5*h);

        for( int i = 0; i < planets.size(); i++){
            K2_r_all.at(i) = planets.at(i).v;
            K2_v_all.at(i) = compute_acceleration(i);
        }
        
        planets = initial_state;

        //k3
        for(int i = 0; i < planets.size(); i++){
            planets.at(i).r += 0.5*h*K2_r_all.at(i); 
            planets.at(i).v += 0.5*h*K2_v_all.at(i);
        }
        //evolveEuler(0.5*h);

        for( int i = 0; i < planets.size(); i++){
            K3_r_all.at(i) = planets.at(i).v;
            K3_v_all.at(i) = compute_acceleration(i);
        }

        planets = initial_state;

        //k4
        for(int i = 0; i < planets.size(); i++){
            planets.at(i).r += h*K3_r_all.at(i); 
            planets.at(i).v += h*K3_v_all.at(i);
        }
        //evolveEuler(h);

        for( int i = 0; i < planets.size(); i++){
            K4_r_all.at(i) = planets.at(i).v;
            K4_v_all.at(i) = compute_acceleration(i);
        }

        for( int i = 0; i < planets.size(); i++){
            final_state.at(i).r = initial_state.at(i).r + (h/6)*(K1_r_all.at(i) + 2*K2_r_all.at(i) + 2*K3_r_all.at(i) + K4_r_all.at(i));
            final_state.at(i).v = initial_state.at(i).v + (h/6)*(K1_v_all.at(i) + 2*K2_v_all.at(i) + 2*K3_v_all.at(i) + K4_v_all.at(i));
        }
        planets = final_state;
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

        //initialize prev_accel as 0
        planets.at(0).a_prev = arma::vec(3).fill(0.);
        planets.at(1).a_prev = arma::vec(3).fill(0.);

    }

double System::adaptive_time_step(double eta, double h){
    //constant eta
    //prev. timestep h
    double min = 10000000000;
    //arbitrarly high number
    for(int j= 0; j<planets.size(); j++){
        double abs_a_over_a_dot;
        double a_dot;
        double a;
        //how to calculate a_dot
        //we need to save the previous pos from every single simulation step 
        a_dot = norm(  (compute_acceleration(j)//System 
                        -planets.at(j).a_prev //System_prev )
                     )/h );
                     //diff. quotient of a 
        
        a = norm(compute_acceleration(j));
        abs_a_over_a_dot = abs(a/a_dot);
        if (eta*abs_a_over_a_dot <min){
            min=eta*abs_a_over_a_dot; 
        }
    } 
    return min;
    //right the timestep becomes massive at eta=3
}