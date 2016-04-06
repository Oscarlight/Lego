/*
 * InOut2D.h
 *
 *  Created on: Nov 12, 2015
 *      Author: oscar
 */

#ifndef SRC_INOUT2D_H_
#define SRC_INOUT2D_H_

#include "GetPot"
#include "Params.h"
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
class InOut2D : public Params {
public:
	/* Read */
	// device general
	std::string devName, userCom;
	int blockNum;
	std::vector<int> leftBT, rightBT;
	std::vector<int> layerNumList;
	// Boundary type
	std::vector<int> topBTList, botBTList;
	// workfunction
	std::vector<double> topWFList, botWFList;

	std::vector<double> blockLength;
	std::vector<int> blockPoint;
	std::vector<int> topGateArea;

	std::vector<std::string> layerName;
	std::vector<double> layerThickness;
	std::vector<int> layerPoint;
	std::vector<bool> connect2Drain;
	std::vector<bool> connect2Source;

	// bias
	double vtgArray[3], vbgArray[3], vdsArray[3], vssArray[3];

	// Poisson Convergence
	double voltageErr;
	double carrierConcenErr; // in cm^-3
	double magicNum;

	// Terminal connection
	bool terminalConnect[3];

	// material
	std::map<std::string, std::vector<double>> matMap;



public:
	InOut2D();
	virtual ~InOut2D();
	// read
	void readDevice2D(std::string filename);
	void readMaterial(std::string filename);
	std::map<int, std::array<double, 4>> gateBiasMap(double Vtg, double Vbg);

	// write
	void writeCB(std::string fileName, std::map<int, std::vector<std::vector<double> > > cBMap);
	void writeBandsinSemi(std::string fileName, std::vector<double> vtgArray, std::vector<double> vdsArray, std::vector<double> vbgArray,
			std::vector< std::vector< std::vector<double> > > band2DPerBias);
	void writeCapaMap(std::string fileName, std::map<std::vector<double>, double> map);
	void writeGateEffMap(std::string fileName, std::map<std::vector<double>, std::vector<double>> map);
	static void printMatToText(const char* filename, mat matrix);
};

#endif /* SRC_INOUT2D_H_ */
