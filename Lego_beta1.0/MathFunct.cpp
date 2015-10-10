/*
 * MathFunct.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: oscar
 */

#include "MathFunct.h"

MathFunct::MathFunct() {
	// TODO Auto-generated constructor stub

}

MathFunct::~MathFunct() {
	// TODO Auto-generated destructor stub
}

double MathFunct::fermiIntegralHalf(double eta)
{
    double result=0;
    double a=0;

    a= pow(eta, 4) + 33.6 * eta * (1 - 0.68 * exp(-0.17* pow((eta +1),2))) + 50;
    result= 2 * sqrt(3.1415926) / ( 3 * sqrt(3.1415926) * pow(a, (-3/8)) + 4* exp(-eta));

    return ( result );
};

double MathFunct::bernoulli(double t)
{
	return ( t/(exp(t)-1) );
};
