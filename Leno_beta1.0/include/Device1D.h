/*
 * MaterialBlock.h
 *
 *  Created on: Oct 8, 2015
 *      Author: oscar
 */

#ifndef DEVICE1D_H_
#define DEVICE1D_H_

// included dependencies
#include "Material.h"

class Device1D : public Material {

protected:
	std::vector<Material> materialList;
	std::vector<int> nyList;
	std::vector<double> dieleArray; // dielectrics in y-direction: size: sumPoint + 1
	std::vector<double> spacingArray; // spacing in y-direction: size: sumPoint + 1
	std::vector<double> electronAffinityArray;
	std::vector<double> bandGapArray;
	std::vector<double> typeArray;
	std::vector<double> phinInitArray; // only service as initial guess
	std::vector<double> phipInitArray;
	int bottomBoundaryType;
	int topBoundaryType;
	int sumPoint;
	sp_mat matrixC;

public:
	Device1D();
	virtual ~Device1D();
	void startWith(Material m, double _t, int _n, int bottomBoundaryType);
	void add(Material m, double _t, int _n); // add 3D material, thickness, how may points inside. If 2D, ignore t, and n
	void endWith(Material m, double _t, int _n, int topBoundaryType);
	void matrixDiff();

	// getter function
	mat getEAArray(); // get electronAffinityArray in mat type
	mat getBGArray(); // get bandGapArray in mat type
	mat getTArray(); // get typeArray in mat type
	mat getPnIArray(); // get phinInitArray in mat type
	mat getPpIArray(); // get phipInitArray in mat type
	std::vector<double> getSpacingY(); // get spacingArray
	int getSumPoint(); // get sumPoint
	int getBBT();
	int getTBT();
	sp_mat getMatrixC();

	// related function
	mat eDensityArrayFunct(mat phin); // mat == Map<double>
	mat hDensityArrayFunct(mat phip);
	mat chargeDensityArrayFunct(mat phin, mat phip, bool Equilibrum);
	sp_mat qCMatFunct(mat phin, mat phip, bool Equilibrum);
	mat vec2Mat(std::vector<double> vec); // change dynamic type vector to static type mat
};

#endif /* DEVICE1D_H_ */
