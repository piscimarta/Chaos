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
}

Planet Planet::operator=(const Planet &other){
    m = other.m;
    r = other.r;
    v = other.v;
    return *this;
}



