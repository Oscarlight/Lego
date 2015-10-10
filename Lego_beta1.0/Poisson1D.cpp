/*
 * Poisson1D.cpp
 *
 *  Created on: Oct 10, 2015
 *      Author: oscar
 */

#include "Poisson1D.h"

Poisson1D::Poisson1D(Device1D _dev1D) {
	dev1D = _dev1D;
	bCArray = zeros(dev1D.sumPoint);
	Vtg = 0; Vbg = 0; WFt = 0; WFb = 0;
	fLnArray = nullptr;
	fLpArray = nullptr;
}

Poisson1D::~Poisson1D() {
	// TODO Auto-generated destructor stub
}

void Poisson1D::setFLnArray(std::vector<double> _fLnArray) {
	fLnArray = _fLnArray;
}

void Poisson1D::setFLpArray(std::vector<double> _fLpArray) {
	fLpArray = _fLpArray;
}

void Poisson1D::setGateBias(double _Vtg, double _WFt, double _Vbg, double _WFb) {
	Vtg = _Vtg; WFt = _WFt; Vbg = _Vbg; WFb = _WFb;
}

mat Poisson1D::getBCArray() {
	return (bCArray);
}
// bCArray initialized in constructor
void Poisson1D::bCAarryFunct() {
    bCArray[0]=(-Vbg + WFb) * LARGE;
    bCArray[dev1D.sumPoint-1]=(-Vtg + WFt) * LARGE;
}
void Poisson1D::run() {
	// TODO how to deal with non-equvilant
}




