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

    // add Constructors  ...    

    Planet();
    Planet operator=(const Planet &other);
};


#endif 
