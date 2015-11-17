/*
 * Device2D.h
 *
 *  Created on: Oct 30, 2015
 *      Author: oscar
 */

#ifndef SRC_DEVICE2D_H_
#define SRC_DEVICE2D_H_

/**
 * Construct the 2D Device based on 1D Device
 * Precondition:
 * 1) uniform spacing in vertical and lateral dicrection
 * 2) each vertical stacks has same material block
 */

#include "Device1D.h"

class Device2D : public Device1D {

private:
	// naming convention: ...List is per 1D block; ...Array is per point
	std::vector<Device1D> dev1DList;
	std::vector<int> nxList;
    std::vector<double> dieleArrayIP; // IP: in-plane, x direction (i.e. lateral)
    std::vector<double> spacingArrayIP;
	std::vector<double> electronAffinityArray;
	std::vector<double> bandGapArray;
	std::vector<double> phinInitArray; // only service as initial guess
	std::vector<double> phipInitArray;
	std::vector<int> leftBTArray;
    std::vector<int> rightBTArray;

	int sumPoint; // totla number of points in device2D
	int unitSize; // total number of vertical points in device1D
	sp_mat matrixC;

public:
	Device2D();
	virtual ~Device2D();
	// Construct Methods
	void startWith_2D(Device1D d, double _w, int _n, std::vector<int> leftBTList);
	void add_2D(Device1D d, double _w, int _n);
	void endWith_2D(Device1D d,  double _w, int _n, std::vector<int> rightBTList);
	void matrixDiff_2D();

	// getter function
	std::vector<Device1D> getDev1DList();
	mat getEAArray_2D(); // get electronAffinityArray in mat type
	mat getBGArray_2D(); // get bandGapArray in mat type
//	mat getTArray_2D(); // get typeArray in mat type
	mat getPnIArray_2D(); // get phinInitArray in mat type
	mat getPpIArray_2D(); // get phipInitArray in mat type
//	std::vector<double> getSpacingY_2D(); // get spacingArray
//	std::vector<double> getSpacingX_2D(); // get spacingArray
	int getSumPoint_2D(); // get sumPoint
	std::vector<int> getNxList(); // getNxList and getUnitSize give all the info of 2D mesh
	int getUnitSize();
//	int getLBT();
//	int getRBT();
	sp_mat getMatrixC_2D();

	// related function
	mat eDensityArrayFunct_2D(mat phin); // mat == Map<double>
	mat hDensityArrayFunct_2D(mat phip);
	mat chargeDensityArrayFunct_2D(mat phin, mat phip, bool Equilibrum);
	sp_mat qCMatFunct_2D(mat phin, mat phip, bool Equilibrum);

// was public because I want to use it also in constructing fLnArray in Poisson2D
// function are not self-contained.
public:
	mat get1DSlice(mat array, int blockIndex, int sliceIndex);
	mat put1DSlice(mat array, mat arraySlice, int blockIndex, int sliceIndex);
	sp_mat put1DSliceInEye(sp_mat eye, sp_mat matSlice, int blockIndex, int sliceIndex);

};

#endif /* SRC_DEVICE2D_H_ */
