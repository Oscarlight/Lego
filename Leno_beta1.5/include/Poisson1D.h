/*
 * Poisson1D.h
 *
 *  Created on: Oct 10, 2015
 *      Author: oscar
 */

#ifndef POISSON1D_H_
#define POISSON1D_H_

// included dependencies
#include "Device1D.h"

class Poisson1D : public Params {
private:
	Device1D dev1D;
	double Vtg;
	double Vbg;
	double WFt;
	double WFb;
	mat fLnArray; // like -Vds array in old version
	mat fLpArray;
	mat bCArray; // boundaryConditionArray
	//
	mat cDA; // chargeDensityArray
	mat mobileED; // mobile eletron density: In transport, sometimes either eletron or hole is involved (contact to either of them)
	mat mobileHD; // mobile hole density
	mat phin;
	mat phip;
	//
	mat condBand;
	mat valeBand;

public:
	Poisson1D(); // default constructor, called when construct its derived class Poisson2D
	Poisson1D(Device1D _dev1D);
	virtual ~Poisson1D();

	mat setFLnArray(std::vector<double> _fLnArray); // if in equilibrium, only need Efn = -V
	mat setFLpArray(std::vector<double> _fLpArray);
	void setGateBias(double _Vtg, double _WFt, double _Vbg, double _WFb);
	void runPoisson1D(double vTolerance, double chargeTolerance, double magicNumber, bool Equilibrum);

	mat getBCArray();
	mat getFLnArray();
	mat getFLpArray();
	mat getCondBand();
	mat getValeBand();
	mat getPhin();
	mat getPhip();
	mat getCDA();
	mat getMobileED(); // - means electron here
	mat getMobileHD(); // - means hole here

	//helper
	std::vector<double> createFLnArray(std::vector<bool> connect2Drain, double Vds,
			std::vector<bool> connect2Source, double Vss);

	std::vector<double> rangeByNum(double begin, double end, double number);
	std::vector<double> rangeByStep(double begin, double end, double step);

protected:
	mat bCArrayFunct(Device1D dev1D, double Vtg, double WFt, double Vbg, double WFb);
	mat setPhip(mat phin, bool Equilibrum);
};

#endif /* POISSON1D_H_ */
