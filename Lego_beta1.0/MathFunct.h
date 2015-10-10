/*
 * MathFunct.h
 *
 *  Created on: Oct 9, 2015
 *      Author: oscar
 */

#ifndef MATHFUNCT_H_
#define MATHFUNCT_H_
#define SQUARE(X) ((X)*(X))

#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <vector>

// include external libs
#include "lib/armadillo-6.100.0/include/armadillo";
#include "lib/cubature-1.0.2/cubature.h";
#include "lib/SuperLU_4.3/SRC/slu_cdefs.h";

using namespace arma;

class MathFunct {
public:
	MathFunct();
	virtual ~MathFunct();
	double fermiIntegralHalf(double eta);
	double bernoulli(double t);
};

#endif /* MATHFUNCT_H_ */
