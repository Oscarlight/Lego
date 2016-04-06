/*
 * MathFunct.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: oscar
 */

#include "MathFunct.h"

MathFunct::MathFunct() {
	// std::cout << "Default constructor of Params is called." << std::endl;

}

MathFunct::~MathFunct() {
	// std::cout << "Default destructor of Params is called." << std::endl;
}

double MathFunct::fermiIntegralHalf(double x)
{
  /* use an analytic expression. The result achieves within 0.4% error in all ranges.*/
  if(x<-4.5)
  {
    //for small arguments, fhfp and exp are almost identical.
    return ( 1.0/exp(-x) );
  }
  else if(x<0.0)
  {
    double v = pow(x,4) + 50 + 33.6*x*(1-0.68*exp(-0.17*(x+1)*(x+1)));
    double p = 1.329340388179*pow(v,double(-0.375));
    return ( 1.0/(exp(-x) + p) );
  }
  else
  {
    double v = pow(x,4) + 50 + 33.6*x*(1-0.68*exp(-0.17*(x+1)*(x+1)));
    double p = 1.329340388179*pow(v,double(-0.375));
    return ( 1.0/(1.0/exp(x) + p) );
  }

};


double MathFunct::bernoulli(double t)
{
	return ( t/(exp(t)-1) );
};
