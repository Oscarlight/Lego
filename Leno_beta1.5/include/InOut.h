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
#include <iostream>
#include <fstream>
#include <sstream>
class InOut : public Params {
public:
	/* Read */
	// device general
	std::string devName, userCom;
	int layerNum;
	// Boundary type
	int topBT, botBT;
	// workfunction
	double topWF, botWF;
	// for plot
	bool plotBA, plotIV, savePl, saveBA, calcuIV;
	double vtgPl[2], vbgPl[2], vdsPl[2];
	// use map don't preserve order
	std::vector<std::string> layerName;
	std::vector<double> layerThickness;
	std::vector<int> layerPoint;
	std::vector<bool> connect2Drain;
	std::vector<bool> connect2Source;
	// bias
	double vtgArray[3], vbgArray[3], vdsArray[3], vssArray[3];
	// material
	std::map<std::string, std::vector<double>> matMap;

	// Poisson Convergence
	double voltageErr;
	double carrierConcenErr; // in cm^-3
	double magicNum;

private:
	std::map<std::vector<double>, std::vector<double>> ivMap;

public:
	InOut();
	virtual ~InOut();

	void readDevice(std::string filename);
	void readMaterial(std::string filename);

	void storeIV(double Vtg, double Vbg, double Vds, std::vector<double> current);
	void writeIV(std::string fileName, std::map<std::vector<double>, std::vector<double>> ivMap);
	void writeBAandCharge(std::string fileName, std::vector<double> vtgArray, std::vector<double> vdsArray, std::vector<double> vbgArray,
			std::vector< std::vector< std::vector<double> > > bandPerBias);

	std::map<std::vector<double>, std::vector<double>> getIVMap();

};

#endif /* INOUT_H_ */
