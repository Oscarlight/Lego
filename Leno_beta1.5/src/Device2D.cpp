/*
 * Device2D.cpp
 *
 *  Created on: Oct 30, 2015
 *      Author: oscar
 */

#include "Device2D.h"

Device2D::Device2D() {
	sumPoint = 0;
	unitSize = 0;
}

Device2D::~Device2D() {
}

void Device2D::startWith_2D(Device1D d, double _w, int _n, std::vector<int> leftBTList) {
	double w, n, spacingX;
	w = _w * cmNL(1E-7);
	n = _n;
	spacingX = w/n;
	/* the total number of device1D points
	 * precondition require this variable to
	 * be the same for each device1D block
	 */
	unitSize = d.getSumPoint();

	/** add left boundary condition to each point in the device1D
	 *
	 */
	for (int i = 0; i < d.getNyList().size(); i++){
		leftBTArray.insert(leftBTArray.end(), d.getNyList()[i], leftBTList[i]);
	}

	dev1DList.push_back(d);
	nxList.push_back(n+1);
	/**
	 * Similar to contruction of 1D device:
	 * 1) add one real point at the left end;
	 * 2) add another fictitious points at the right end.
	 */
	// TODO: It was + 2, but here I change to +1
	for (int i = 0; i < n + 1; i++){
		/* dieleArrayIP.insert(dieleArrayIP.end(), d.getDieleArrayIP().begin(), d.getDieleArrayIP().end());
		 * 	Wrong !! : d.getDieleArrayIP() return "different" copy of the same vector!
		*/
		std::vector<double> tempV0 = d.getDieleArrayIP();
		// std::cout << "DieleArray " << tempV0.size() << std::endl;
		dieleArrayIP.insert(dieleArrayIP.end(), tempV0.begin(), tempV0.end());

		spacingArrayIP.insert(spacingArrayIP.end(), unitSize, spacingX);

		std::vector<double> tempV1 = d.getPnIArrayVec();
		// std::cout << "PniArray " << tempV1.size() << std::endl;
		phinInitArray.insert(phinInitArray.end(), tempV1.begin(), tempV1.end());

		std::vector<double> tempV2 = d.getPpIArrayVec();
		phipInitArray.insert(phipInitArray.end(), tempV2.begin(), tempV2.end());

		std::vector<double> tempV3 = d.getBGArrayVec();
		bandGapArray.insert(bandGapArray.end(), tempV3.begin(), tempV3.end());

		std::vector<double> tempV4 = d.getEAArrayVec();
		electronAffinityArray.insert(electronAffinityArray.end(), tempV4.begin(), tempV4.end());
	}
	std::cout << "1D Device Block is added with # of lateral point = " << n + 1
			<< ", and vertical point = " << unitSize << std::endl;
}

void Device2D::add_2D(Device1D d, double _w, int _n) {
	double w, n, spacingX;
	w = _w * cmNL(1E-7);
	n = _n;
	spacingX = w/n;

	if (d.getSumPoint() != unitSize)
		std::cerr << "Error: Each 1D Block Should have the same vertical points!" << std::endl;

	dev1DList.push_back(d);
	nxList.push_back(n);

	//TODO: here was n + 2, follow the original code but I forgot why
	// here I use n
	for (int i = 0; i < n; i++){
		std::vector<double> tempV0 = d.getDieleArrayIP();
		dieleArrayIP.insert(dieleArrayIP.end(), tempV0.begin(), tempV0.end());

		spacingArrayIP.insert(spacingArrayIP.end(), unitSize, spacingX);

		std::vector<double> tempV1 = d.getPnIArrayVec();
		phinInitArray.insert(phinInitArray.end(), tempV1.begin(), tempV1.end());

		std::vector<double> tempV2 = d.getPpIArrayVec();
		phipInitArray.insert(phipInitArray.end(), tempV2.begin(), tempV2.end());

		std::vector<double> tempV3 = d.getBGArrayVec();
		bandGapArray.insert(bandGapArray.end(), tempV3.begin(), tempV3.end());

		std::vector<double> tempV4 = d.getEAArrayVec();
		electronAffinityArray.insert(electronAffinityArray.end(), tempV4.begin(), tempV4.end());
	}
	std::cout << "1D Device Block is added with # of lateral point = " << n
			<< ", and vertical point = " << unitSize << std::endl;
}

void Device2D::endWith_2D(Device1D d,  double _w, int _n, std::vector<int> rightBTList) {
	/**
	 * add right boundary conditions to all layers
	 */
	for (int i = 0; i < d.getNyList().size(); i++){
		rightBTArray.insert(rightBTArray.end(), d.getNyList()[i], rightBTList[i]);
	}
	add_2D(d, _w, _n);
	// total # of read points
	sumPoint = dieleArrayIP.size();
	// add another one fictitious points at the right end

	std::vector<double> tempV = d.getDieleArrayIP();
	dieleArrayIP.insert(dieleArrayIP.end(), tempV.begin(), tempV.end());
	spacingArrayIP.insert(spacingArrayIP.end(), unitSize, spacingArrayIP.back());
	std::cout << "Device2D is constructed with total points = " << sumPoint << std::endl;
}

void Device2D::matrixDiff_2D() {
	int numBlock = dev1DList.size();
	// initialize the matrixC2D;
	mat matrixC2D(dev1DList[0].getMatrixC()); // has to make it dense mat to use join_rows/join_cols
	// TODO: was j= 0
	for (int j = 0; j < numBlock; j++) {
		std::cout << "Constructing differential matrix @ Block #" << j+1 << "; sumPoint = " << dev1DList[j].getSumPoint()
				<< " * " << nxList[j] << std::endl;
		dev1DList[j].matrixDiff();
		mat matrixC1D(dev1DList[j].getMatrixC());
		// std::cout << matrixC1D.size() << std::endl;
		for (int i=0; i < nxList[j]; i++) {
			matrixC2D=join_cols( join_rows(matrixC2D, zeros(matrixC2D.n_rows, dev1DList[j].getSumPoint())),
					join_rows(zeros(dev1DList[j].getSumPoint(), matrixC2D.n_cols), matrixC1D) );
		}
	}

	// initialize the C_L
	// was sumPointLateral * unitSize, here sumPoint = sumPointLateral * unitSize
	mat C_L = zeros(sumPoint, sumPoint);
	// similar to Device1D 1 => unitSize
    for ( int i=unitSize; i< sumPoint + unitSize; i++) {
        for ( int j=unitSize; j< sumPoint + unitSize; j++) {
            if (i==j) {
                C_L(i-unitSize,j-unitSize)= - 2*dieleArrayIP[i]/( spacingArrayIP[i]*( spacingArrayIP[i]+spacingArrayIP[i-unitSize] ) )
                 - 2*dieleArrayIP[i-unitSize]/( spacingArrayIP[i]*( spacingArrayIP[i]+spacingArrayIP[i-unitSize] ) );
            } else if (j==i-unitSize) {
                C_L(i-unitSize,j-unitSize)=2*dieleArrayIP[i-unitSize]/( spacingArrayIP[i]*( spacingArrayIP[i]+spacingArrayIP[i-unitSize] ) );
            } else if (j==i+unitSize) {
                C_L(i-unitSize,j-unitSize)=2*dieleArrayIP[i]/( spacingArrayIP[i]*( spacingArrayIP[i]+spacingArrayIP[i-unitSize] ) );
            }
        }
    }
    // Lateral boundary
    // boundary condition can differ at the semi and at the dielectric, apply point by point at right and left edge.

    /* IMPORTANT: Ohmic -> Neumann; Schottchy -> Dirichlet
     *
     * TODO: implement Schottchy: 1. add workfucntion for contact metal in input file;
     * 							  2. add this workfunction to boundary;
     * 							  3. Vds apply to boundary condition array instead of VdsArray (see notebook IV)
    */

    for (int i = 0; i < unitSize; i++) {
    	// left boundary
    	if ( leftBTArray[i] == Dirichlet ) {
    		C_L(i,unitSize + i)=LARGE;
    	} else if ( leftBTArray[i] == Neumann) {
    		C_L(i,unitSize + i)= - C_L(i,i);
    	} else {
    		std::cerr << "Error: left boundary type not found." << std::endl;
    	}
    	// right boundary
    	if ( rightBTArray[i] == Dirichlet ) {
    		C_L(i + sumPoint - unitSize, i + sumPoint - 2*unitSize) = LARGE;
    	} else if ( rightBTArray[i] == Neumann) {
    		C_L(i + sumPoint - unitSize, i + sumPoint - 2*unitSize) = -C_L(i + sumPoint - unitSize, i + sumPoint - unitSize);
    	} else {
    		std::cerr << "Error: right boundary type not found." << std::endl;
    	}
    }
    // std::cout << matrixC2D.size() << ", " << C_L.size() << std::endl;
    sp_mat tempMatrixC(matrixC2D + C_L);
    matrixC = tempMatrixC;
}

std::vector<Device1D> Device2D::getDev1DList() {
	return ( dev1DList );
}

mat Device2D::getEAArray_2D() {
	return ( electronAffinityArray );
}

mat Device2D::getBGArray_2D() {
	return ( bandGapArray );
}

//mat Device2D::getTArray_2D() {
//}
//
mat Device2D::getPnIArray_2D() {
	return ( phinInitArray );
}

mat Device2D::getPpIArray_2D() {
	return ( phipInitArray );
}
//
//std::vector<double> Device2D::getSpacingY_2D() {
//}
//
std::vector<double> Device2D::getSpacingX_2D() {
	return (spacingArrayIP);
}

std::vector<int> Device2D::getNxList() {
	return ( nxList );
}

int Device2D::getUnitSize() {
	return ( unitSize );
}

int Device2D::getSumPoint_2D() {
	return ( sumPoint );
}

//int Device2D::getLBT() {
//}
//
//int Device2D::getRBT() {
//}

sp_mat Device2D::getMatrixC_2D() {
	return ( matrixC );
}


mat Device2D::chargeDensityArrayFunct_2D(mat phin, mat phip, bool Equilibrum) {
	mat density = zeros(sumPoint);
	for (int i = 0; i < nxList.size(); i++) { // for each 2D device block
		for (int j = 0; j < nxList[i]; j++) { // for each 1D slice in each 2D block
			mat densitySlice = dev1DList[i].chargeDensityArrayFunct(get1DSlice(phin, i, j), get1DSlice(phip, i, j), Equilibrum);
			// std::cout << "densitySlice size: " << densitySlice.size() << std::endl;
			density = put1DSlice(density, densitySlice, i, j);
		}
	}
	return (density);
}

sp_mat Device2D::qCMatFunct_2D(mat phin, mat phip, bool Equilibrum) {
	sp_mat qc(sumPoint, sumPoint); // bug 11/16: speye -> zeros
	for (int i = 0; i < nxList.size(); i++) { // for each 2D device block
		for (int j = 0; j < nxList[i]; j++) { // for each 1D slice in each 2D block
			sp_mat qcSlice =  dev1DList[i].qCMatFunct(get1DSlice(phin, i, j), get1DSlice(phip, i, j), Equilibrum);
			// std::cout << "qcSlice: " << qcSlice << std::endl;
			// std::cout << "qcSlice size: " << qcSlice.size() << std::endl;
			qc = put1DSliceInEye(qc, qcSlice, i, j);
		}
	}
	return ( qc );
}

mat Device2D::get1DSlice(mat array, int blockIndex, int sliceIndex) {
	if (array.n_elem != sumPoint)
		std::cerr << "Error: input arrays of put1DSlice are of wrong size. Expected: " <<
		sumPoint << ". Get: " << array.n_elem << std::endl;

	mat a = zeros(unitSize);
	int sumSliceNum = 0;
	for (int i = 0; i < blockIndex; i++) {
		sumSliceNum += nxList[i];
	}
	int accumIndex = (sumSliceNum + sliceIndex)*unitSize;
	for (int j = 0; j < unitSize; j++) {
		a(j) = array(accumIndex + j);
	}
	return ( a );
}

mat Device2D::put1DSlice(mat array, mat arraySlice, int blockIndex, int sliceIndex) {
	if (array.n_elem != sumPoint || arraySlice.n_elem != unitSize)
		std::cerr << "Error: input arrays of put1DSlice are of wrong size. Expected: (" <<
		sumPoint << ", " << unitSize << "). Get: (" <<
		array.n_elem << ", " << arraySlice.n_elem << ")" << std::endl;

	int sumSliceNum = 0;
	for (int i = 0; i < blockIndex; i++) {
		sumSliceNum += nxList[i];
	}
	int accumIndex = (sumSliceNum + sliceIndex)*unitSize;
	for (int j = 0; j < unitSize; j++) {
		array(accumIndex + j) = arraySlice(j); // reserved assignment compared to get1DSlice
	}
	return ( array );
}


sp_mat Device2D::put1DSliceInEye(sp_mat eye, sp_mat matSlice, int blockIndex, int sliceIndex) {
	if (eye.size() != sumPoint*sumPoint || matSlice.n_elem != unitSize * unitSize)
		std::cerr << "Error: input eye matrix or slice array are of wrong size. Expected: (" <<
		sumPoint*sumPoint << ", " << unitSize * unitSize << "). Get: (" <<
		eye.size() << ", " << matSlice.n_elem << ")" << std::endl;

	int sumSliceNum = 0;
	for (int i = 0; i < blockIndex; i++) {
		sumSliceNum += nxList[i];
	}
	int accumIndex = (sumSliceNum + sliceIndex)*unitSize;
	for (int j = 0; j < unitSize; j++) {
		eye(accumIndex + j, accumIndex + j) = matSlice(j,j); // reserved assignment compared to get1DSlice // bug 11/16 (j) -> (j)(j);
	}
	return ( eye );
}
