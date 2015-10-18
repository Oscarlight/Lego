/*
 * InOut.h
 *
 *  Created on: Oct 15, 2015
 *      Author: oscar
 */

#ifndef INOUT_H_
#define INOUT_H_

#include "GetPot"
#include "Params.h"
#include <map>
class InOut : public Params {
public:
	// device general
	std::string devName, userCom;
	int layerNum;
	// Boundary type
	int topBT, botBT;
	// workfunction
	double topWF, botWF;
	// for plot
	bool plotBA, plotIV, saveBA, savePl;
	double vtgPl[2], vbgPl[2], vdsPl[2];
	// use map don't preserve order
	std::vector<std::string> layerName;
	std::vector<double> layerThickness;
	std::vector<int> layerPoint;

	double vtgArray[3], vbgArray[3], vdsArray[3];
	// material
	std::map<std::string, std::vector<double>> matMap;
public:
	InOut();
	virtual ~InOut();

	void readDevice(std::string filename);
	void readMaterial(std::string filename);
};

#endif /* INOUT_H_ */
