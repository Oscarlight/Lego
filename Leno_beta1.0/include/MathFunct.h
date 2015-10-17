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
/**
 * SuperLU for sparse matrix
 * uncomment #define ARMA_USE_SUPERLU in config.hpp in /include/armadillo-bits
 */
#define ARMA_USE_SUPERLU // come before #include <armadillo>
#include "armadillo"
#include "slu_cdefs.h"
#include "cubature.h"

using namespace arma;

class MathFunct {
public:
	MathFunct();
	virtual ~MathFunct();
	double fermiIntegralHalf(double eta);
	double bernoulli(double t);
};

#endif /* MATHFUNCT_H_ */
