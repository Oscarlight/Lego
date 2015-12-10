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
	// for 1D
	std::vector<double> x, spacing, cB, vB, fLn, fLp, chargeDensity;
	mat phin, phip;
	std::vector<double> xSemi, cBSemi, vBSemi, fLnSemi, fLpSemi, mcSemi, mvSemi;
	// for 2D
	/**
	 * in the map, key is the block index; first vector index is the slice index, second vector index is the point index in a slice
	 *  (slice: 1D device)
	 */
	std::map<int, std::vector<std::vector<double> > > cBMap, vBMap, chargeDensityMap, spacingXMap;
	std::vector<double> topGateChargeArray; // 1/cm^2
private:
	mat typeArray;
	bool bandAndChargeDone;
	bool bandAndCharge2DDone;
public:
	ExtractData();;
	virtual ~ExtractData();

	bool bandAndCharge(Poisson1D p1D, Device1D dev1D);
	bool bandAndCharge2D(Poisson2D p2D, Device2D dev2D);

	double topGateCharge2D(Device2D dev2D, std::vector<int> gateArea);

	/* band alignment in semiconductor only
	 * alway run after bandAlignment
	 */
	int bASemiOnly();

	std::vector<double> mat2Vec(mat m);
	std::map<int, std::vector<std::vector<double>>> array2Map(mat a, int unitSize, std::vector<int> nxList);
	std::vector<double> mapRowSlide(std::map<int, std::vector<std::vector<double> > > map, int rowIndex);
	double arrayPoint(std::vector<double> array, int pointIndex);
};

#endif /* EXTRACTDATA_H_ */
