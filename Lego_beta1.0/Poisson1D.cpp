/*
 * Poisson1D.cpp
 *
 *  Created on: Oct 10, 2015
 *      Author: oscar
 */

#include "Poisson1D.h"

Poisson1D::Poisson1D(Device1D _dev1D) {
	dev1D = _dev1D;
	bCArray = zeros(dev1D.getSumPoint());
	Vtg = 0; Vbg = 0; WFt = 0; WFb = 0;
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

void Poisson1D::bCArrayFunct() {
	// bCArray initialized in constructor
	// if boundary type is Neumann, all set to zero
	if (dev1D.getBBT() == Dirichlet)
		bCArray[0]=(-Vbg + WFb) * LARGE;
	if (dev1D.getTBT() == Dirichlet)
		bCArray[dev1D.getSumPoint()-1]=(-Vtg + WFt) * LARGE;
}

mat Poisson1D::setPhip(mat phin, bool Equilibrum) {
	if (Equilibrum == true)
		return ( zeros(dev1D.getSumPoint()) ); // see Device1D::chargeDensityArrayFunct

	return ( dev1D.getBGArray() - phin - (dev1D.vec2Mat(fLnArray) - dev1D.vec2Mat(fLpArray)) );
}

void Poisson1D::runPoisson1D(double vTolerance, double _chargeTolerance, double magicNumber, bool Equilibrum) {
	// Relatively quicker to solve potential, because potential is continuous.
	// potential == vacuum level
	double chargeTolerance = _chargeTolerance / ( pow(cmNL(1), 3) );
	mat fLnArrayM = dev1D.vec2Mat(fLnArray);
	mat potential= phin + dev1D.getEAArray() + fLnArrayM;
	bCArrayFunct();
	mat oldcDA = zeros(dev1D.getSumPoint());
	mat errV = ones(1)*LARGE;
	mat errCDA = ones(1)*LARGE;
	// superlu_opts settings; // optional
    int numConvergenceStep=0;
    do {
    	cDA = dev1D.chargeDensityArrayFunct(phin, setPhip(phin, Equilibrum), Equilibrum);
    	errCDA = abs(cDA - oldcDA).max();
    	oldcDA = cDA;
    	mat error = - (dev1D.getMatrixC() * potential - ChargeQ/E0 * cDA  - bCArray);
    	sp_mat qCMat = ChargeQ/E0 * dev1D.qCMatFunct(phin, setPhip(phin, Equilibrum), Equilibrum);
    	sp_mat matrixC_plusCq = dev1D.getMatrixC() - qCMat;
    	mat deltaPotential = spsolve(matrixC_plusCq, error, "superlu");
    	/** spsolve was not recognized. solve by using full path in #include,
    	 *  then follow http://stackoverflow.com/questions/30494610/how-to-link-armadillo-with-eclipse
    	 *  then rebuild (in index)
    	 */
    	potential += magicNumber * deltaPotential;
    	errV = abs(deltaPotential).max();
    	numConvergenceStep++;
    	phin = potential - (dev1D.getEAArray() + fLnArrayM );
    } while ((errV[0] > vTolerance) || (errCDA[0] > chargeTolerance));
    phip = setPhip(phin, Equilibrum);
}




