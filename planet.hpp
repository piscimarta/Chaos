// include guard
#ifndef __planet_hpp__   
#define __planet_hpp__

// include headers
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <armadillo>


class Planet
{
public:
    // variables
    double m;
    arma::vec r = arma::vec(3);
    arma::vec v = arma::vec(3);
    arma::vec a_prev = arma::vec(3);
    //save a planets previous acceleration in the planet class

    // add Constructors  ...    

    Planet();
    Planet operator=(const Planet &other);
};


#endif 
