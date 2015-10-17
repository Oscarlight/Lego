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
protected:
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
	mat phin;
	mat phip;
	//
	mat condBand;
	mat valeBand;

public:
	Poisson1D(Device1D _dev1D);
	virtual ~Poisson1D();

	void setFLnArray(std::vector<double> _fLnArray); // if in equilibrium, only need Efn = -V
	void setFLpArray(std::vector<double> _fLpArray);
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

	std::vector<double> rangeByNum(double begin, double end, double number);
	std::vector<double> rangeByStep(double begin, double end, double step);

protected:
	void bCArrayFunct();
	mat setPhip(mat phin, bool Equilibrum);
};

#endif /* POISSON1D_H_ */
