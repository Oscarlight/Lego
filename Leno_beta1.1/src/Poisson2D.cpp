/*
 * Poisson2D.cpp
 *
 *  Created on: Nov 4, 2015
 *      Author: oscar
 */

#include "Poisson2D.h"

Poisson2D::Poisson2D(Device2D _dev2D) {
	dev2D = _dev2D;
	bCArray = zeros(dev2D.getSumPoint_2D());
	phin = dev2D.getPnIArray_2D();
	phip = dev2D.getPnIArray_2D();
}

Poisson2D::~Poisson2D() {
}

void Poisson2D::setFLnArray_2D(std::map<int, std::vector<double>> fLnMap) {
	fLnArray = zeros(dev2D.getSumPoint_2D());
	// validating
	if (dev2D.getNxList().size() != fLnMap.size())
		std::cerr << "Error: electron Fevel Level Map failed to cover all 1D Blocks!" << std::endl;

	for (int i = 0; i < dev2D.getNxList().size(); i++) { // index of 1D block
		for (int j = 0; j < dev2D.getNxList()[i]; j++) { // index of 1D slice in a 1d block // bug 11/16: dev2D.getNxList()[i]: j -> i
//			std::cout << "fLnArray: " << setFLnArray(fLnMap.at(i)) << std::endl;
			fLnArray = dev2D.put1DSlice(fLnArray, setFLnArray(fLnMap.at(i)), i, j); // bug 11/16: forget to assign back
		}
	}
}

void Poisson2D::setFLpArray_2D(std::map<int, std::vector<double>> fLpMap) {
	fLpArray = zeros(dev2D.getSumPoint_2D());
	// validating
	if (dev2D.getNxList().size() != fLpMap.size())
		std::cerr << "Error: hole Fevel Level Map failed to cover all 1D Blocks!" << std::endl;

	for (int i = 0; i < dev2D.getNxList().size(); i++) { // index of 1D block
		for (int j = 0; j < dev2D.getNxList()[i]; j++) { // index of 1D slice in a 1d block // bug 11/16: dev2D.getNxList()[i]: j -> i
			fLpArray = dev2D.put1DSlice(fLpArray, setFLpArray(fLpMap.at(i)), i, j);  // bug 11/16: forget to assign back
		}
	}
}

void Poisson2D::setGateBias_2D(std::map<int, std::array<double, 4>> _biasMap) {
	// validating
	if (dev2D.getNxList().size() != _biasMap.size())
		std::cerr << "Error: bias map failed to cover all 1D Blocks!" << std::endl;

	biasMap = _biasMap;
}

mat Poisson2D::bCArrayFunct_2D(Device2D dev2D, std::map<int, std::array<double, 4>> biasMap) {
	mat bCArray = zeros(dev2D.getSumPoint_2D());

	for (int i = 0; i < dev2D.getNxList().size(); i++) { // index of 1D block
		for (int j = 0; j < dev2D.getNxList()[i]; j++) { // index of 1D slice in a 1d block // bug 11/16: dev2D.getNxList()[i]: j -> i
			double tempVtg = biasMap.at(i)[0];
			double tempWFt = biasMap.at(i)[1];
			double tempVbg = biasMap.at(i)[2];
			double tempWFb = biasMap.at(i)[3];
//			std::cout << "bCArray: " << bCArrayFunct(dev2D.getDev1DList()[i], tempVtg, tempWFt, tempVbg, tempWFb) << std::endl;
			bCArray = dev2D.put1DSlice(bCArray, bCArrayFunct(dev2D.getDev1DList()[i],
					tempVtg, tempWFt, tempVbg, tempWFb), i, j);
		}
	}
	this->bCArray = bCArray;
	return ( bCArray );
}

mat Poisson2D::getBCArray_2D() {
	return (bCArray);
}

mat Poisson2D::getFLnArray_2D() {
	return (fLnArray);
}

mat Poisson2D::getFLpArray_2D() {
	return (fLpArray);
}

mat Poisson2D::getCondBand_2D() {
	return (condBand);
}

mat Poisson2D::getValeBand_2D() {
	return (valeBand);
}

mat Poisson2D::getPhin_2D() {
	return (phin);
}

mat Poisson2D::getPhip_2D() {
	return (phip);
}

mat Poisson2D::getCDA_2D() {
	return (cDA);
}

mat Poisson2D::setPhip_2D(mat phin, bool Equilibrum) {
	if (Equilibrum == true)
		return ( zeros(dev2D.getSumPoint_2D()) );

	return ( dev2D.getBGArray_2D() - phin - (fLnArray - fLpArray));
}

void Poisson2D::runPoisson2D(double vTolerance, double _chargeTolerance, double magicNumber, bool Equilibrum) {
	// Relatively quicker to solve potential, because potential is continuous.
	// potential == vacuum level
	double chargeTolerance = _chargeTolerance / ( pow(cmNL(1), 3) );
	// std::cout << phin.size() << ", " << dev2D.getEAArray_2D().size() << ", " << fLnArray.size() << std::endl;
	mat potential= phin + dev2D.getEAArray_2D() + fLnArray;
//  std::cout << dev2D.getEAArray_2D() << std::endl; // proved!
//	std::cout << fLnArray << std::endl; /* Mark 11/16: found bugs in fLnArray */
	mat bCArray = bCArrayFunct_2D(dev2D, biasMap);
//	std::cout << bCArray << std::endl; /* Mark 11/16: found bugs in bCArray */
	mat oldcDA = zeros(dev2D.getSumPoint_2D());
	mat errV = ones(1)*LARGE;
	mat errCDA = ones(1)*LARGE;
	// superlu_opts settings; // optional
	// mat denseMat(dev2D.getMatrixC_2D());
	// InOut2D::printMatToText("matrixC", denseMat);
    int numConvergenceStep=0;
    do {
    	cDA = dev2D.chargeDensityArrayFunct_2D(phin, setPhip_2D(phin, Equilibrum), Equilibrum);
    	errCDA = abs(cDA - oldcDA).max();
    	oldcDA = cDA;
//    	 std::cout << cDA << std::endl;
    	// std::cout << dev2D.getMatrixC_2D().size() << ", " << potential.size() << ", " << cDA.size() << ", " << bCArray.size() << std::endl;
    	mat error = - (dev2D.getMatrixC_2D() * potential - ChargeQ/E0 * cDA  - bCArray);
//    	 std::cout << error << std::endl;
    	sp_mat qCMat = ChargeQ/E0 * dev2D.qCMatFunct_2D(phin, setPhip_2D(phin, Equilibrum), Equilibrum);
//    	 std::cout << qCMat << std::endl;
    	sp_mat matrixC_plusCq = dev2D.getMatrixC_2D() - qCMat; // TODO: all use 2D version, use Private field help avoid errors!
    	mat deltaPotential = spsolve(matrixC_plusCq, error, "superlu");

    	/** spsolve was not recognized. solve by using full path in #include,
    	 *  then follow http://stackoverflow.com/questions/30494610/how-to-link-armadillo-with-eclipse
    	 *  then rebuild (in index)
    	 */

//    	std::cout << phin << std::endl;

		potential += magicNumber * deltaPotential;
    	errV = abs(deltaPotential).max();
    	numConvergenceStep++;
    	phin = potential - (dev2D.getEAArray_2D() + fLnArray );

    	// std::cout << errV(0) << std::endl;

    } while ((errV(0) > vTolerance) || (errCDA(0) > chargeTolerance) || numConvergenceStep < 1E4);
    if ( numConvergenceStep > 1E4)
    	std::cerr << "+++++ Solution not found: convergence step over 1E4 +++++" << std::endl;
    phip = setPhip_2D(phin, Equilibrum);
    condBand = potential - dev2D.getEAArray_2D();
    valeBand = condBand - dev2D.getBGArray_2D();
}
