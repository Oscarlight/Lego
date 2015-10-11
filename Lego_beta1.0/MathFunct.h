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
// #include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

// include external libs
// #include "/home/oscar/Documents/CPP/workspace/Lego_beta1.0/lib/armadillo-6.100.0/include/armadillo";
// #include "/home/oscar/Documents/CPP/workspace/Lego_beta1.0/lib/cubature-1.0.2/cubature.h";
/**
 * SuperLU for sparse matrix
 * uncomment #define ARMA_USE_SUPERLU in config.hpp in /include/armadillo-bits
 */
// #include "/home/oscar/Documents/CPP/workspace/Lego_beta1.0/lib/SuperLU_4.3/SRC/slu_cdefs.h";

/*
 *
 */
#include <armadillo>
#define ARMA_USE_SUPERLU
#include "lib/SuperLU_4.3/SRC/slu_cdefs.h"
#include "lib/cubature-1.0.2/cubature.h"

using namespace arma;

class MathFunct {
public:
	MathFunct();
	virtual ~MathFunct();
	double fermiIntegralHalf(double eta);
	double bernoulli(double t);
};

#endif /* MATHFUNCT_H_ */
