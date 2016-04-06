/*
 * Poisson1D.cpp
 *
 *  Created on: Oct 10, 2015
 *      Author: oscar
 */

#include "Poisson1D.h"

Poisson1D::Poisson1D() {
}

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

std::vector<double> Poisson1D::createFLnArray(std::vector<bool> connect2Drain,
		double Vds, std::vector<bool> connect2Source, double Vss) {
	std::vector<double> fL(dev1D.getSumPoint(), 0); // initalize to 0

	int interAccu = 0;
	for (int i = 0; i < dev1D.getNyList().size(); i++) { // for each layer
		// std::cout << dev1D.getNyList()[i] << std::endl;
		// std::cout << connect2Drain[i] << "; " << connect2Source[i] << std::endl;
		if (connect2Drain[i] == true) {
			for (int j = 0; j < dev1D.getNyList()[i]; j++) { // for each point in a layer
				fL[interAccu + j] = -Vds;
			}
		}
		if (connect2Source[i] == true) {
			for (int j = 0; j < dev1D.getNyList()[i]; j++) { // for each point in a layer
				fL[interAccu + j] = -Vss;
			}
		}
		interAccu += dev1D.getNyList()[i];
	}

	return (fL);
}


mat Poisson1D::setFLnArray(std::vector<double> _fLnArray) {
	fLnArray = dev1D.vec2Mat(_fLnArray);
	return ( fLnArray ); // add return for Poisson2D
}

mat Poisson1D::setFLpArray(std::vector<double> _fLpArray) {
	fLpArray = dev1D.vec2Mat(_fLpArray);
	return ( fLpArray );
}

void Poisson1D::setGateBias(double _Vtg, double _WFt, double _Vbg, double _WFb) {
	Vtg = _Vtg; WFt = _WFt; Vbg = _Vbg; WFb = _WFb;
}

mat Poisson1D::getBCArray() {
	return (bCArray);
}

mat Poisson1D::bCArrayFunct(Device1D dev1D, double Vtg, double WFt, double Vbg, double WFb) {
	// bCArray initialized in constructor
	// if boundary type is Neumann, all set to zero
	// keep variable as self-constain as possible
	mat locBCArray = zeros(dev1D.getSumPoint());
	if (dev1D.getBBT() == Dirichlet)
		locBCArray(0)=(-Vbg + WFb) * LARGE;
	if (dev1D.getTBT() == Dirichlet)
		locBCArray(dev1D.getSumPoint()-1)=(-Vtg + WFt) * LARGE;
	bCArray = locBCArray; // TODO: if bCArray has not been used before bCArrayFunct, we can initialize it here instead of in constructor
	return ( bCArray );
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

mat Poisson1D::getMobileED() {
	return (mobileED);
}

mat Poisson1D::setPhip(mat phin, bool Equilibrum) {
	if (Equilibrum == true)
		return (dev1D.getBGArray() - phin);
		// see Device1D::chargeDensityArrayFunct

	return ( dev1D.getBGArray() - phin - (fLnArray - fLpArray) );
}

void Poisson1D::runPoisson1D(double vTolerance, double _chargeTolerance, double magicNumber, bool Equilibrum) {
	// Relatively quicker to solve potential, because potential is continuous.
	// potential == vacuum level
	double chargeTolerance = _chargeTolerance / ( pow(cmNL(1), 3) );
	mat potential= phin + dev1D.getEAArray() + fLnArray;
	bCArrayFunct(dev1D, Vtg, WFt, Vbg, WFb);
	mat oldcDA = zeros(dev1D.getSumPoint());
	mat errV = ones(1)*LARGE;
	mat errCDA = ones(1)*LARGE;
	// superlu_opts settings; // optional
    int numConvergenceStep=0;
//    std::cout << " error has NaN: " << setPhip(phin, Equilibrum) << std::endl;
    do {
    	cDA = dev1D.chargeDensityArrayFunct(phin, setPhip(phin, Equilibrum), Equilibrum);
    	errCDA = abs(cDA - oldcDA).max();
    	oldcDA = cDA;

    	// std::cout << cDA(51) << " "  << cDA(87) << std::endl;
    	mat error = - (dev1D.getMatrixC() * potential - ChargeQ/E0 * cDA  - bCArray);

    	sp_mat qCMat = ChargeQ/E0 * dev1D.qCMatFunct(phin, setPhip(phin, Equilibrum), Equilibrum);
    	sp_mat matrixC_plusCq = dev1D.getMatrixC() - qCMat;

    	/** TODO: debug the matrix
    	 */
    	// mat matrixC_plusCq_dense(matrixC_plusCq);
    	// mat deltaPotential = solve(matrixC_plusCq_dense, error);

    	mat deltaPotential = spsolve(matrixC_plusCq, error, "superlu");
        /** spsolve(): could not solve system
         *  	t2D in material file must be equal to the spacingY!
         *  	for now, if set to t2D = 0.1, magicNumber = 1, endless loop
         */

    	/** spsolve was not recognized. solve by using full path in #include,
    	 *  then follow http://stackoverflow.com/questions/30494610/how-to-link-armadillo-with-eclipse
    	 *  then rebuild (in index)
    	 */
    	potential += magicNumber * deltaPotential;
    	errV = abs(deltaPotential).max();
    	numConvergenceStep++;
    	phin = potential - (dev1D.getEAArray() + fLnArray );
    } while ((errV(0) > vTolerance) || (errCDA(0) > chargeTolerance) || numConvergenceStep < 1E4);
    if ( numConvergenceStep > 1E5)
    	std::cerr << "+++++ Solution not found: convergence step over 1E5 +++++" << std::endl;
    phip = setPhip(phin, Equilibrum);
    condBand = potential - dev1D.getEAArray();
    valeBand = condBand - dev1D.getBGArray();
    mobileED =  - dev1D.mobileEDensityArrayFunct(phin); // positive means hole
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
	do {
		v.push_back(begin);
		begin += step;
	} while(begin * 1E5 < ceil((end + step)*1E5));
	// v.push_back(end);
	return (v);
}


