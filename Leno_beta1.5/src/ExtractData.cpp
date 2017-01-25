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
	potentialMap = array2Map(p2D.getPotential_2D(), dev2D.getUnitSize(), dev2D.getNxList());
	cBMap = array2Map(p2D.getCondBand_2D(), dev2D.getUnitSize(), dev2D.getNxList());
	vBMap = array2Map(p2D.getValeBand_2D(), dev2D.getUnitSize(), dev2D.getNxList());
	chargeDensityMap = array2Map(p2D.getCDA_2D(), dev2D.getUnitSize(), dev2D.getNxList()); //still in Leno units
	spacingXMap = array2Map(dev2D.getSpacingX_2D(), dev2D.getUnitSize(), dev2D.getNxList());
	return (bandAndCharge2DDone = true);
}


/** return C/cm not charge density */
double ExtractData::topGateCharge2D(Device2D dev2D, std::vector<int> gateArea) {
	sumCharge = 0;
	// std::cout << " Enter top gate charge 2D" << std::endl;
	for (int i : gateArea) {
		// std::cout << " gateArea = " << i <<  std::endl;
		std::vector<std::vector<double>> cB, spacingX;
		cB = cBMap.at(i);
		spacingX = spacingXMap.at(i);
		double dieleY = dev2D.getDev1DList()[i].getMaterialList().back().dieleY;
		double gateThick = dev2D.getDev1DList()[i].getNyList().back() * dev2D.getDev1DList()[i].getSpacingY().back() * cmLN(1); // in cm
		double gateCap = dieleY * ( E0 * cLN(1)/cmLN(1)) / gateThick;
		// charge density (1/cm^2) at each points
		for (int j = 0; j < cB.size(); j++) {
			// std::cout << cB[j].size() - dev2D.getDev1DList()[i].getNyList().back() << ", " << cB[j].size() << std::endl;
			// std::cout << cB[j][ cB[j].size() - dev2D.getDev1DList()[i].getNyList().back() ] << ", " << cB[j].back() << std::endl;
			double vOverTopOx = cB[j][ cB[j].size() - dev2D.getDev1DList()[i].getNyList().back() ] - cB[j].back();
			// std::cout << "vOverTopOx = " << vOverTopOx << std::endl;
			double tempCharge = gateCap * vOverTopOx; // C/cm^2
			// std::cout << gateCap << ", " << vOverTopOx << ", " << tempCharge << std::endl;
			topGateChargeArray.push_back(tempCharge);
			sumCharge += tempCharge * ( spacingX[j].back() * cmLN(1) ); // C/cm
		}
	}
	return ( sumCharge );
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
	mobileElectronDensity = mat2Vec(p1D.getMobileED() / ( cmLN(1) * cmLN(1) * cmLN(1) ) );
	spacing = dev1D.getSpacingY();
	double sum = 0;
	for (int i = 0; i < cB.size(); i++) {
		x.push_back(sum);
		sum += (cmLN(spacing[i])*1E7); // x in nm;
	}
	typeArray = dev1D.getTArray();
	return (bandAndChargeDone = true);
}

/**
 * Run it after bandAndCharge
 */
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

/**
 * Precondition: Run it after bandAndCharge2D
 * write to bandByLayer and chargeByLayer,
 * which have band alignment (conduction band/valence band)
 * and charge density in 2D layers.
 * return band by layer for the layers connected to S/D
 */
// TODO: generalize it, now design only for pin/thin tfet
int ExtractData::bAChargeSemiOnly2D(std::vector<int> sdIndex){
	if( bandAndCharge2DDone != true)
		std::cerr << "Error: bandAlignment 2D needs to be run first." << std::endl;

	int numOfBand = 2; // for pin-TFET, only 2: conduction and valence
	bandByLayer.push_back(mapRowSlide(cBMap, sdIndex[0])); // top for ThinTFET
	bandByLayer.push_back(mapRowSlide(vBMap, sdIndex[0]));
	chargeByLayer.push_back(mapRowSlide(chargeDensityMap, sdIndex[0]));
	if (sdIndex[0] != sdIndex[1]) { // source and drain are not at the same layer
		bandByLayer.push_back(mapRowSlide(cBMap, sdIndex[1])); // bottom for ThinTFET
		bandByLayer.push_back(mapRowSlide(vBMap, sdIndex[1]));
		chargeByLayer.push_back(mapRowSlide(chargeDensityMap, sdIndex[1]));
		numOfBand += 2; // for Thin-TFET, 2 more bands
	}
	return (numOfBand);
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
	return ( map );
}

// row slice, i.e. slice in x-direction
std::vector<double> ExtractData::mapRowSlide(std::map<int, std::vector<std::vector<double> > > map, int rowIndex) {
	std::vector<double> row;
	int accu = 0;
	for (int i = 0; i < map.size(); i++) { // per block
		std::vector<std::vector<double> > block = map.at(i);
		for (int j = 0; j < block.size(); j++) { // per slice in block
			std::vector<double> slice = block[j];
			row.push_back(slice[rowIndex]);
		}
	}
	return ( row );
}

double ExtractData::arrayPoint(std::vector<double> array, int pointIndex) {
	return ( array[pointIndex] );
}
