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
	std::map<int, std::vector<std::vector<double> > > cBMap, vBMap, chargeDensityMap;
private:
	mat typeArray;
	bool bandAndChargeDone;
	bool bandAndCharge2DDone;
public:
	ExtractData();;
	virtual ~ExtractData();

	bool bandAndCharge(Poisson1D p1D, Device1D dev1D);
	bool bandAndCharge2D(Poisson2D p2D, Device2D dev2D);
	/* band alignment in semiconductor only
	 * alway run after bandAlignment
	 */
	int bASemiOnly();

	std::vector<double> mat2Vec(mat m);
	std::map<int, std::vector<std::vector<double>>> array2Map(mat a, int unitSize, std::vector<int> nxList);

};

#endif /* EXTRACTDATA_H_ */
