
// include guard
#ifndef __header_hpp__
#define __header_hpp__

// include headers
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>

// Return a string with a double in scientific notation
std::string scientific_format(double d, const int &width, const int &prec);

// Return a string with a vector<double> in scientific notation
std::string scientific_format(const std::vector<double> &v, const int &width, const int &prec);

//Return a string of an integer
std::string int_to_str(int d);


#endif // end of include guard 
