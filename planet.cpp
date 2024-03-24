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

// include header
#include "planet.hpp"



Planet::Planet()
{
    // assign values of charge, mass, position and velocity
    m = 1;
    r = arma::vec(3).fill(0.);
    v = arma::vec(3).fill(0.);
    a_prev = arma::vec(3).fill(0.);
    //save the previous acc. in the planet class 
}

Planet Planet::operator=(const Planet &other){
    m = other.m;
    r = other.r;
    v = other.v;
    a_prev = other.a_prev;
    return *this;
}


// here we need to define the constructor of the class

