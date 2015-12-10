/*
 * Capacitance.h
 *
 *  Created on: Nov 20, 2015
 *      Author: oscar
 */

#ifndef SRC_CAPACITANCE_H_
#define SRC_CAPACITANCE_H_

#include "Params.h"
#include <map>

class Capacitance : public Params {

private:
	std::vector<double> vtgArray;
	std::vector<double> vbgArray;
	std::vector<double> vdsArray;
	std::map<std::vector<double>, double> topGateChargeMap, cgsMap, cgdMap, cggMap;

public:
	Capacitance();
	virtual ~Capacitance();

	std::map<double, double> cgsMeyer();
	std::map<double, double> cgdMeyer();

	void readCharge(double Vtg, double Vbg, double Vds, double Vss, double charge);
	void calCggCgsCgd(std::vector<double> vtgArray, std::vector<double> vbgArray, std::vector<double> vssArray, std::vector<double> vdsArray, bool terminalConnect[3]);

	std::map<std::vector<double>, double> getCgsMap();
	std::map<std::vector<double>, double> getCgdMap();
	std::map<std::vector<double>, double> getCggMap();
};

#endif /* SRC_CAPACITANCE_H_ */
