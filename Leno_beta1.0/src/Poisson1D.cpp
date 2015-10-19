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
	phin = dev1D.getPnIArray();
	phip = dev1D.getPnIArray();
	// std::cout << "Default constructor of Poisson1D is called." << std::endl;
}

Poisson1D::~Poisson1D() {
	// std::cout << "Default destructor of Poisson1D is called." << std::endl;
}

void Poisson1D::setFLnArray(std::vector<double> _fLnArray) {
	fLnArray = dev1D.vec2Mat(_fLnArray);
}

void Poisson1D::setFLpArray(std::vector<double> _fLpArray) {
	fLpArray = dev1D.vec2Mat(_fLpArray);
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

mat Poisson1D::getFLnArray() {
	return ( fLnArray );
}

mat Poisson1D::getFLpArray() {
	return ( fLpArray );
}

mat Poisson1D::getCondBand() {
	return ( condBand );
}

mat Poisson1D::getValeBand() {
	return ( valeBand );
}

mat Poisson1D::getPhin() {
	return (phin);
}

mat Poisson1D::getPhip() {
	return (phip);
}

mat Poisson1D::getCDA() {
	return (cDA);
}

mat Poisson1D::setPhip(mat phin, bool Equilibrum) {
	if (Equilibrum == true)
		return ( zeros(dev1D.getSumPoint()) ); // see Device1D::chargeDensityArrayFunct

	return ( dev1D.getBGArray() - phin - (fLnArray - fLpArray) );
}

void Poisson1D::runPoisson1D(double vTolerance, double _chargeTolerance, double magicNumber, bool Equilibrum) {
	// Relatively quicker to solve potential, because potential is continuous.
	// potential == vacuum level
	double chargeTolerance = _chargeTolerance / ( pow(cmNL(1), 3) );
	mat potential= phin + dev1D.getEAArray() + fLnArray;
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
    	phin = potential - (dev1D.getEAArray() + fLnArray );
    } while ((errV(0) > vTolerance) || (errCDA(0) > chargeTolerance));
    phip = setPhip(phin, Equilibrum);
    condBand = potential - dev1D.getEAArray();
    valeBand = condBand - dev1D.getBGArray();
}

std::vector<double> Poisson1D::rangeByNum(double begin, double end, double number) {
	std::vector<double> v;
	double it = (end - begin)/(number - 1);
	v.push_back(begin);
	for (int i = 0; i < number; i++){
		begin += it;
		v.push_back(it);
	}
	return (v);
}

std::vector<double> Poisson1D::rangeByStep(double begin, double end, double step) {
	std::vector<double> v;
	// Preventing double precision error (< 10000)
	for (double i = begin; i*10000 < ceil((end)*10000); i+=step){
		// std::cout << i*10000 << ceil((end-step)*10000) << std::endl;
		v.push_back(i);
	}
	v.push_back(end);
	return (v);
}
