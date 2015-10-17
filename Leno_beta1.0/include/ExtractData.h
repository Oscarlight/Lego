/*
 * Plot.h
 *
 *  Created on: Oct 11, 2015
 *      Author: oscar
 */

#ifndef EXTRACTDATA_H_
#define EXTRACTDATA_H_


#include "Poisson1D.h"

class ExtractData : public Params {

public:
	std::vector<double> x, spacing, cB, vB, fLn, fLp, chargeDensity;
	mat phin, phip;
	std::vector<double> xSemi, cBSemi, vBSemi, fLnSemi, fLpSemi, mcSemi, mvSemi;
private:
	mat typeArray;
	bool bandAndChargeDone;
public:
	ExtractData();;
	virtual ~ExtractData();

	bool bandAndCharge(Poisson1D p1D, Device1D dev1D);

	/* band alignment in semiconductor only
	 * alway run after bandAlignment
	 */
	int bASemiOnly();

	std::vector<double> mat2Vec(mat m);
};

#endif /* EXTRACTDATA_H_ */
