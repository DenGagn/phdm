/*!
* \file GrapheneFloquet.cpp
*
*
* \author Author: D. Gagnon <denisg6@hotmail.com>
*/

// Include some headers
#include <iostream>
#include <fstream>
#include <armadillo>
#include <cmath>
#include <complex>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "GammaFactor.hpp"
#include "Utils.hpp"

/// Main function
int main(int argc, char *argv[])
{

    auto Gamma = GammaValues(2.30,0.0,316227766016837.94,1e10,2);
    
    for(auto it = Gamma.begin(); it != Gamma.end(); it++)
    {
        std::cout << *it << std::endl;
    }
    
}
