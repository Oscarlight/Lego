/*
 * Plot.h
 *
 *  Created on: Oct 11, 2015
 *      Author: oscar
 */

#ifndef EXTRACTDATA_H_
#define EXTRACTDATA_H_


#include "Poisson1D.h"
#include "Poisson2D.h"

class ExtractData : public Params {

public:
	// ---- for 1D band information
	std::vector<double> x, spacing, cB, vB, fLn, fLp, chargeDensity, mobileElectronDensity;
	mat phin, phip;
	std::vector<double> xSemi, cBSemi, vBSemi, fLnSemi, fLpSemi, mcSemi, mvSemi;

	// ---- for 2D band information
	/**
	 * in the map, key is the block index; first vector index is the slice index, second vector index is the point index in a slice
	 *  (slice: 1D device)
	 */
	std::map<int, std::vector<std::vector<double> > > cBMap, vBMap, chargeDensityMap, spacingXMap;
	/**
	 * organized in conduction band 1, valence band 1, conduction band 2, valence band 2, ...
	 */
	std::vector<std::vector<double> > bandByLayer;

	// ---- for capacitance
	// charge density array on the top gate, set in topGateCharge2D()
	std::vector<double> topGateChargeArray; // 1/cm^2
	// total charge on the top gate, also set in topGateCharge2D()
	double sumCharge; // C/cm

private:
	mat typeArray;
	bool bandAndChargeDone;
	bool bandAndCharge2DDone;

public:
	ExtractData();;
	virtual ~ExtractData();

	// ---- for 1D band info
	bool bandAndCharge(Poisson1D p1D, Device1D dev1D);
	// ---- for 2D band info
	bool bandAndCharge2D(Poisson2D p2D, Device2D dev2D);
	// ---- for capacitance
	double topGateCharge2D(Device2D dev2D, std::vector<int> gateArea);

	// ---- for 1D tunnelling calculation
	/* band alignment in semiconductor only
	 * alway run after bandAlignment
	 * return: num of Semi
	 */
	int bASemiOnly();

	// ---- for 2D capacitance calculation
	/* band alignment in selected layers
	 * TODO: find a better way to generize it
	 */
	int bASemiOnly2D(std::vector<int> sdIndex);


	// Helper methods
	std::vector<double> mat2Vec(mat m);
	std::map<int, std::vector<std::vector<double>>> array2Map(mat a, int unitSize, std::vector<int> nxList);
	std::vector<double> mapRowSlide(std::map<int, std::vector<std::vector<double> > > map, int rowIndex);
	double arrayPoint(std::vector<double> array, int pointIndex);
};

#endif /* EXTRACTDATA_H_ */
