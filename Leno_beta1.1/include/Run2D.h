/*
 * Run2D.h
 *
 *  Created on: Nov 5, 2015
 *      Author: oscar
 */

#ifndef SRC_RUN2D_H_
#define SRC_RUN2D_H_

#include <stdio.h>
#include "Device2D.h"
#include "ExtractData.h"
#include "Poisson2D.h"
#include "InOut2D.h"


class Run2D {
private:
	InOut2D io2D;
	std::vector<ExtractData> band2DPerBias;

public:
	Run2D();
	virtual ~Run2D();
	std::map<std::string, Material> readInput(int argc, char** argv);
	Device2D createDevice2D(int argc, char** argv);
	void runPoisson2D(int argc, char** argv);

	std::vector<ExtractData> getBand2D();
	InOut2D getIO2D();

private:
	// i: index of the block, accu: accumalted number of layers in previous blocks
	Device1D createDevice1D(InOut2D io2D, std::map<std::string, Material> matLib, int i, int accu);
};

#endif /* SRC_RUN2D_H_ */
