#include <pybind11/pybind11.h>
#include <boost/version.hpp>
#include <iostream>
#include <stdlib.h>

namespace py=pybind11;

void print_boost_version(){
std::cout << "Using Boost "     
          << BOOST_VERSION / 100000     << "."  // major version
          << BOOST_VERSION / 100 % 1000 << "."  // minor version
          << BOOST_VERSION % 100                // patch level
          << std::endl;

}


PYBIND11_MODULE(tocheck_boostversion, m) {
    
    m.def("print_boost_version", &print_boost_version);
    
}
