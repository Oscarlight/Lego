/*
 * Plot.cpp
 *
 *  Created on: Oct 11, 2015
 *      Author: oscar
 */

#include "ExtractData.h"

ExtractData::ExtractData() {
	// std::cout << "Default constructor of Plot is called." << std::endl;
}

ExtractData::~ExtractData() {
	// std::cout << "Default destructor of Plot is called." << std::endl;
}

std::vector<double> ExtractData::mat2Vec(mat m) {
	std::vector<double> vec;
	for (int i = 0; i < m.n_elem; i ++) {
		vec.push_back(m(i));
	}
	return (vec);
}

bool ExtractData::bandAndCharge(Poisson1D p1D, Device1D dev1D) {
	// std::vector<double> x, spacing, cB, vB, fLn, fLp;
	cB = mat2Vec(p1D.getCondBand());
	vB = mat2Vec(p1D.getValeBand());
	fLn = mat2Vec(p1D.getFLnArray());
	fLp = mat2Vec(p1D.getFLpArray());
	phin = p1D.getPhin();
	phip = p1D.getPhip();
	chargeDensity = mat2Vec(p1D.getCDA() / ( cmLN(1) * cmLN(1) * cmLN(1) ) );
	spacing = dev1D.getSpacingY();
	double sum = 0;
	for (int i = 0; i < cB.size(); i++) {
		x.push_back(sum);
		sum += cmLN(spacing[i]);
	}
	typeArray = dev1D.getTArray();
	return (bandAndChargeDone = true);
}


int ExtractData::bASemiOnly() {
	if( bandAndChargeDone != true)
		std::cerr << "Error: bandAlignment needs to be run first." << std::endl;
	int numSemi = 0;
	for(int i = 0; i < typeArray.n_elem; i++) {
		if( typeArray(i) == Semiconductor ) {
			xSemi.push_back(x[i]);
			cBSemi.push_back(cB[i]);
			vBSemi.push_back(vB[i]);
			fLnSemi.push_back(fLn[i]);
			fLpSemi.push_back(fLp[i]);
			mcSemi.push_back(phin(i));
			mvSemi.push_back(phip(i));
			numSemi ++;
		}
	}
	return (numSemi);
}
