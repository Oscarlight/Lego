/*
 * run1D.h
 *
 *  Created on: Nov 4, 2015
 *      Author: oscar
 */

#ifndef SRC_RUN1D_H_
#define SRC_RUN1D_H_

/**
 * run Poisson 1D
 */
#include <stdio.h>
#include "Device1D.h"
#include "ExtractData.h"
#include "Poisson1D.h"
#include "Tunnelling.h"
#include "InOut.h"
using namespace std;
class Run1D {
private:
	InOut io;
	std::vector<ExtractData> bandPerBias; // store the band diagram info per bias

public:
	Run1D();
	virtual ~Run1D();
	std::map<std::string, Material> readInput(int argc, char** argv);
	Device1D createDevice1D(int argc, char** argv);
	void runPoisson1D(int argc, char** argv);

	//getter
	InOut getIO();
	std::vector<ExtractData> getBand();
};

#endif /* SRC_RUN1D_H_ */
