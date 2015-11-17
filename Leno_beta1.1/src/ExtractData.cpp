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

bool ExtractData::bandAndCharge2D(Poisson2D p2D, Device2D dev2D) {
	cBMap = array2Map(p2D.getCondBand_2D(), dev2D.getUnitSize(), dev2D.getNxList());
	vBMap = array2Map(p2D.getValeBand_2D(), dev2D.getUnitSize(), dev2D.getNxList());
	chargeDensityMap = array2Map(p2D.getCDA_2D(), dev2D.getUnitSize(), dev2D.getNxList());
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
		sum += (cmLN(spacing[i])*1E7); // x in nm;
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

std::map<int, std::vector<std::vector<double> > > ExtractData::array2Map(mat a, int unitSize, std::vector<int> nxList) {
	int accu = 0;
	std::map<int, std::vector<std::vector<double> > > map;
	for (int i = 0; i < nxList.size(); i++) { // per block
		std::vector<std::vector<double> > block;
		for (int j = 0; j < nxList[i]; j++) { // per slice in block
			std::vector<double> slice;
			for (int k = 0; k < unitSize; k++) { // per point in slice
				slice.push_back(a(accu*unitSize + k));
			}
			accu++;
			block.push_back(slice);
		}
		map.insert(std::pair<int, std::vector<std::vector<double> >>(i, block));
	}
	return map;
}
