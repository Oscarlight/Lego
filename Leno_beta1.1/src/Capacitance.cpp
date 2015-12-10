/*
 * Capacitance.cpp
 *
 *  Created on: Nov 20, 2015
 *      Author: oscar
 */

#include "Capacitance.h"

Capacitance::Capacitance() {
	// TODO Auto-generated constructor stub

}

Capacitance::~Capacitance() {
	// TODO Auto-generated destructor stub
}

void Capacitance::readCharge(double Vtg, double Vbg, double Vds, double Vss, double charge){
	std::vector<double> v;
	v.push_back(Vtg), v.push_back(Vbg), v.push_back(Vds), v.push_back(Vss);
	// std::cout << "Charge = " << charge << std::endl;
	topGateChargeMap.insert(std::pair<std::vector<double>, double>(v, charge));
}

void Capacitance::calCggCgsCgd(std::vector<double> vtgArray, std::vector<double> vbgArray, std::vector<double> vdsArray, std::vector<double> vssArray, bool terminalConnect[3]) {
	for (int k = 0; k < vbgArray.size(); k++) {
		double Vbg = vbgArray[k];

		// Cgg for all Vds, Vtg
		/* for Cgg, Vss should always = 0; */

		for (int i = 0; i < vdsArray.size(); i++) {
			double Vds = vdsArray[i];
			// std::cout << "Calculate Cgd at Vds = " << Vds << std::endl;
			if ( vtgArray.size() > 1 ) {
				for (int j = 0; j < vtgArray.size()-1; j++) {
					double nextVtg = vtgArray[j+1];
					double thisVtg = vtgArray[j];

					if (terminalConnect[0] == true) // Vtg == Vbg
						Vbg = nextVtg;

					std::vector<double> nextV; nextV.push_back(nextVtg); nextV.push_back(Vbg); nextV.push_back(Vds); nextV.push_back(vssArray[0]);

					if (terminalConnect[0] == true) // Vtg == Vbg
						Vbg = thisVtg;

					std::vector<double> thisV; thisV.push_back(thisVtg); thisV.push_back(Vbg); thisV.push_back(Vds); thisV.push_back(vssArray[0]);
					double nextCharge = topGateChargeMap.at(nextV);
					double thisCharge = topGateChargeMap.at(thisV);
					double dvtg = nextVtg - thisVtg;
					double dCharge = nextCharge - thisCharge;
					double Cgs = dCharge / dvtg;
					// std::cout << "--- Cgs = " << Cgs << " @ Vtg = " << thisVtg << std::endl;
					// Vss == 0, so exclude it from voltage array for print (i.e. thisVPrint)
					std::vector<double> thisVPrint; thisVPrint.push_back(Vbg); thisVPrint.push_back(Vds); thisVPrint.push_back(thisVtg); // change order
					cggMap.insert(std::pair<std::vector<double>, double>(thisVPrint, Cgs));
				}
			}
		}
		std::cout << "***** Cgg map is generated. " << std::endl;


		// Cgd for all Vtg, Vds
		/* for Cgd, Vss should always = 0; */

		for (int i = 0; i < vtgArray.size(); i++) {
			double Vtg = vtgArray[i]; // Vtg == Vgs
			// std::cout << "Calculate Cgd at Vtg = " << Vtg << std::endl;
			if ( vdsArray.size() > 1 ) {
				for (int j = 0; j < vdsArray.size()-1; j++) {
					// std::cout << "Calculate Cgd at Vds = " << vdsArray[j] << std::endl;
					// std::cout << "Calculate Cgd at Next Vds = " << vdsArray[j+1] << std::endl;
					double nextVgd = Vtg - vdsArray[j+1]; // Vtg - nextVds = Vgs - Vds = Vgd;
					double thisVgd = Vtg - vdsArray[j]; // Vtg - thisVds;

					if (terminalConnect[0] == true) // Vtg == Vbg
						Vbg = Vtg;

					std::vector<double> nextV; nextV.push_back(Vtg); nextV.push_back(Vbg); nextV.push_back(vdsArray[j+1]); nextV.push_back(vssArray[0]);
					std::vector<double> thisV; thisV.push_back(Vtg); thisV.push_back(Vbg); thisV.push_back(vdsArray[j]); thisV.push_back(vssArray[0]);
					double nextCharge = topGateChargeMap.at(nextV);
					double thisCharge = topGateChargeMap.at(thisV);
					double dvgd = nextVgd - thisVgd;
					double dCharge = nextCharge - thisCharge;
					double Cgd = dCharge / dvgd;
					// std::cout << "--- Cgd = " << Cgd << " @ Vgd = " << thisVgd << std::endl;
					// Vss == 0, so exclude it from voltage array for print (i.e. thisVPrint)
					std::vector<double> thisVPrint; thisVPrint.push_back(Vbg); thisVPrint.push_back(Vtg); thisVPrint.push_back(vdsArray[j]);

					if (terminalConnect[0] == true) // if double gated, dOg = 2 dOtg
						Cgd = 2 * Cgd;

					cgdMap.insert(std::pair<std::vector<double>, double>(thisVPrint, Cgd));
				}
			}
		}
		std::cout << "***** Cgd map is generated. " << std::endl;

		// Cgs for all Vtg, Vds
		/* for Cgs, Vss will be two values, just used to cal Cgs; */
		if ( vssArray.size() != 2 || vssArray[0] != 0.0)
			std::cerr << "Vss Array Size is not 2, or the frist Vss is not 0.0." << std::endl;
		else {
			for (int i = 0; i < vtgArray.size(); i++) {
				double Vtg = vtgArray[i]; // Vtg == Vgs
				for (int j = 0; j < vdsArray.size(); j++) {

					if (terminalConnect[0] == true) // Vtg == Vbg
						Vbg = Vtg;

					double Vds = vdsArray[j];

					std::vector<double> nextV; nextV.push_back(Vtg); nextV.push_back(Vbg); nextV.push_back(Vds); nextV.push_back(vssArray[1]);
					std::vector<double> thisV; thisV.push_back(Vtg); thisV.push_back(Vbg); thisV.push_back(Vds); thisV.push_back(vssArray[0]);
					double nextCharge = topGateChargeMap.at(nextV);
					double thisCharge = topGateChargeMap.at(thisV);
					double dvss = vssArray[1] - vssArray[0];
					double dCharge = nextCharge - thisCharge;
					double Cgs = dCharge / dvss;
					// std::cout << "--- Cgd = " << Cgd << " @ Vgd = " << thisVgd << std::endl;
					// Vss == 0, so exclude it from voltage array for print (i.e. thisVPrint)
					std::vector<double> thisVPrint; thisVPrint.push_back(Vbg); thisVPrint.push_back(Vtg); thisVPrint.push_back(Vds);
					cgsMap.insert(std::pair<std::vector<double>, double>(thisVPrint, Cgs));
				}
			}
			std::cout << "***** Cgs map is generated. " << std::endl;
		}


		if (terminalConnect[0] == true) // Vtg == Vbg
			break;
	}

}

std::map<std::vector<double>, double> Capacitance::getCgsMap() {
	return ( cgsMap );
}

std::map<std::vector<double>, double> Capacitance::getCgdMap() {
	return ( cgdMap );
}

std::map<std::vector<double>, double> Capacitance::getCggMap() {
	return ( cggMap );
}
