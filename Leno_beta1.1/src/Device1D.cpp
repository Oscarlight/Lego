/*
 * MaterialBlock.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: oscar
 */

#include "Device1D.h"

Device1D::Device1D() {
	sumPoint = 0;
	topBoundaryType = 0;
	bottomBoundaryType = 0;
	// std::cout << "Default constructor of Device1D is called." << std::endl;
}

Device1D::~Device1D() {
	// std::cout << "Default destructor of Device1D is called." << std::endl;
}

/**
 * start building the device1D
 * 	initialize the following arrays: using vector<double> (like ArrayList in Java)
 * 		for storing all the materials:
 * 			materialList
 * 		for MatrixC:
 * 			dieleArray
 * 			spacingArray
 * 		for electron Affinity Array:
 * 			electronAffinityArray
 * 		for bandGap Array:
 * 			bandGapArray
 * 		for material type Array:
 * 			typeArray
 * 		for phin and phip initialization:
 * 			phinArray
 * 			phipArray
 *
 */
void Device1D::startWith(Material m, double _t, int _n, int _bottomBoundaryType) {
	double t, n;
	double spacingY;

	// 3D or 2D
	if (m.dimension == 3) {
		t = _t * cmNL(1E-7);
		n = _n;
		spacingY = t/n;
	} else if ( m.dimension == 2) {
		t = m.t2D;
		n = 1;
		spacingY = t;
	}

	// bottomBonudaryType
	bottomBoundaryType = _bottomBoundaryType;

	// store the material and the # of read points
	materialList.push_back(m);
	nyList.push_back(n+1);
	typeList.push_back(m.type);
	/**
    * take an example
    * Original: Ny: 3, 3, 2, 3, 3: # of points = 14
    *  |     SiO2  |bottomSemi | vdWGap | topSemi   |    SiO2   |
    *      *   *   *   #   #   #   @    @   #   #   #   *   *   *
    * Add one real point at the bottom end: so dieleArray.size = # of points + 1 = 15 = sumPoint
    *  |     SiO2  |bottomSemi | vdWGap | topSemi   |    SiO2   |
    *  *   *   *   *   #   #   #   @    @   #   #   #   *   *   *
    * In order to constuct the finite difference method, add another one fictitious points at the end: 16
    *  |     SiO2  |bottomSemi | vdWGap | topSemi   |    SiO2   |
    *  *   *   *   *   #   #   #   @    @   #   #   #   *   *   *   *
    * e[0]e[1]e[2]e[3]...                                          e[15]
    * dieleArray is descreatized in the middle of two points, so set e(i-1/2) -> e(i-1); e(i+1/2) -> e(i)
	*/

	// build arrays
	// add one read point at the bottom end
	dieleArray.insert(dieleArray.end(), n + 1, m.dieleY);
	dieleArrayIP.insert(dieleArrayIP.end(), n + 1, m.dieleX); // for 2D later

	spacingArray.insert(spacingArray.end(), n + 1, spacingY);
	electronAffinityArray.insert(electronAffinityArray.end(), n + 1, m.electronAffinity);
	bandGapArray.insert(bandGapArray.end(), n + 1, m.bandGap);
	typeArray.insert(typeArray.end(), n + 1, m.type);
	( m.type == Semiconductor ) ? phinInitArray.insert(phinInitArray.end(), n + 1, m.bandGap/2) : phinInitArray.insert(phinInitArray.end(), n + 1, 0);
	( m.type == Semiconductor ) ? phipInitArray.insert(phipInitArray.end(), n + 1, m.bandGap/2) : phipInitArray.insert(phipInitArray.end(), n + 1, 0);

	std::cout << "Material " << m.name << " is added as first layer with # of point = " << n + 1 << std::endl;
}

void Device1D::add(Material m, double _t, int _n) {
	double t, n;
	double spacingY;

	// 3D or 2D
	// m.dimension is prohibited, if dimension is protected.
	if (m.dimension == 3) {
		t = _t * cmNL(1E-7);
		n = _n;
		spacingY = t/n;
	} else if (m.dimension == 2) {
		t = m.t2D;
		n = 1;
		spacingY = t;
	}
	// store the material and the # of read points
	materialList.push_back(m);
	nyList.push_back(n);
	typeList.push_back(m.type);
	// 	extend arrays
	dieleArray.insert(dieleArray.end(), n, m.dieleY);
	dieleArrayIP.insert(dieleArrayIP.end(), n, m.dieleX); // for 2D later

	spacingArray.insert(spacingArray.end(), n, spacingY);
	electronAffinityArray.insert(electronAffinityArray.end(), n, m.electronAffinity);
	bandGapArray.insert(bandGapArray.end(), n, m.bandGap);
	typeArray.insert(typeArray.end(), n, m.type);
	( m.type == Semiconductor ) ? phinInitArray.insert(phinInitArray.end(), n, m.bandGap/2) : phinInitArray.insert(phinInitArray.end(), n, 0);
	( m.type == Semiconductor ) ? phipInitArray.insert(phipInitArray.end(), n, m.bandGap/2) : phipInitArray.insert(phipInitArray.end(), n, 0);

	std::cout << "Material " << m.name << " is added with # of point = " << n << std::endl;
}

void Device1D::endWith(Material m, double _t, int _n, int _topBoundaryType) {
	// topBoundaryType
	topBoundaryType = _topBoundaryType;
	add(m, _t, _n);
	// add another one fictitious points at the end
	dieleArray.push_back(m.dieleY);
	// dieleArrayIP don't need to add one!

	spacingArray.push_back(spacingArray.back());
	// total # of read points
	sumPoint = dieleArray.size() - 1;
	std::cout << " *** Device1D is constructed with total points = " << sumPoint << std::endl;
}



mat Device1D::vec2Mat(std::vector<double> vec) {
	mat t = zeros(vec.size());
	for (int i = 0; i < vec.size(); i++) {
		t(i) = vec[i];  // ALWAYS use () for mat type
	}
	return ( t );
}

void Device1D::matrixDiff() {
	sp_mat C(sumPoint, sumPoint);
	mat dieleArrayM = vec2Mat(dieleArray);
	mat spacingArrayM = vec2Mat(spacingArray);
	/**
	 * Diff Matrix considering B.C. using Penelty Method
	 */
	for ( int i=1; i<sumPoint+1; i++) {
		for ( int j=1; j<sumPoint+1; j++) {
			if (i==j) {
				C(i-1,j-1)= - 2*dieleArrayM(i)/( spacingArrayM(i)*( spacingArrayM(i)+spacingArrayM(i-1) ) )
						- 2*dieleArrayM(i-1)/( spacingArrayM(i)*( spacingArrayM(i)+spacingArrayM(i-1) ) );
			} else if (j==i-1) {
				C(i-1,j-1)=2*dieleArrayM(i-1)/( spacingArrayM(i)*( spacingArrayM(i)+spacingArrayM(i-1) ) );
			} else if (j==i+1) {
				C(i-1,j-1)=2*dieleArrayM(i)/( spacingArrayM(i)*( spacingArrayM(i)+spacingArrayM(i-1) ) );
			};
		};
	};
	// Here both end can have different boundaries type
	if ( bottomBoundaryType == Dirichlet ) {
		C(0,0)=LARGE;   // panelty method to enforse boundary condition
	} else if (bottomBoundaryType == Neumann) {
		C(0,0)= - C(0,1);
	} else {
		std::cerr << "Error: bottom boundary type not found." << std::endl;
	}

	if ( topBoundaryType == Dirichlet ) {
		C(sumPoint-1,sumPoint-1)=LARGE; // panelty method to enforse boundary condition
	} else if (topBoundaryType == Neumann) {
		C(sumPoint-1,sumPoint-1)= - C(sumPoint-1,sumPoint-2);
	} else {
		std::cerr << "Error: top boundary type not found." << std::endl;
	}
	matrixC = C;
}


std::vector<Material> Device1D::getMaterialList() {
	return ( materialList );
}

mat Device1D::getEAArray() {
	return ( vec2Mat(electronAffinityArray) );
}
mat Device1D::getBGArray() {
	return ( vec2Mat(bandGapArray) );
}
std::vector<double> Device1D::getEAArrayVec() {
	return ( electronAffinityArray );
}
std::vector<double> Device1D::getBGArrayVec() {
	return ( bandGapArray );
}
mat Device1D::getTArray() {
	return ( vec2Mat(typeArray) );
}
mat Device1D::getPnIArray() {
	return ( vec2Mat(phinInitArray) );
}
mat Device1D::getPpIArray() {
	return ( vec2Mat(phipInitArray) );
}
std::vector<double> Device1D::getPnIArrayVec() {
	return ( phinInitArray );
}
std::vector<double> Device1D::getPpIArrayVec() {
	return ( phipInitArray );
}
int Device1D::getSumPoint() {
	return ( sumPoint );
}
int Device1D::getBBT() {
	return ( bottomBoundaryType);
}
int Device1D::getTBT() {
	return ( topBoundaryType);
}
sp_mat Device1D::getMatrixC() {
	return ( matrixC );
}
std::vector<double> Device1D::getSpacingY() {
	return ( spacingArray );
}
std::vector<double> Device1D::getDieleArrayIP() {
	return ( dieleArrayIP );
}
std::vector<int> Device1D::getNyList() {
	return ( nyList );
}

mat Device1D::eDensityArrayFunct(mat phin) {
	if ( phin.n_elem != sumPoint)
		std::cerr << "Error: input wrong phin size!" << std::endl;
    mat density=zeros(sumPoint);

    int k=0;
    for ( int i=0; i < materialList.size(); i++) {
        for ( int j=0; j < nyList[i]; j++) {
            if (typeList[i] == Semiconductor) {
                density(k) = materialList[i].electronDensity( (double)phin(k) ); // why here used to be k+1, then "carrierDensity[0]=carrierDensity[1];"
            } else if (typeList[i]==Dielectric) {
             };
            k++;
        };
    };
    return ( density );
}

mat Device1D::hDensityArrayFunct(mat phip) {
	if ( phip.n_elem != sumPoint)
		std::cerr << "Error: input wrong phip size!" << std::endl;
    mat density=zeros(sumPoint);

    int k=0;
    for ( int i=0; i < materialList.size(); i++) {
        for ( int j=0; j < nyList[i]; j++) {
            if (typeList[i] == Semiconductor) {
                density(k) = materialList[i].holeDensity( (double)phip(k) ); // why here used to be k+1, then "carrierDensity[0]=carrierDensity[1];"
            } else if (typeList[i]==Dielectric) {
             };
            k++;
        };
    };
    return ( density );
}

mat Device1D::chargeDensityArrayFunct(mat phin, mat phip, bool Equilibrum) {
	if ( phin.n_elem != sumPoint || ( phip.n_elem != sumPoint ) )
		std::cerr << "Error: input wrong phin/phip size!" << std::endl;
    mat density=zeros(sumPoint);

	if (Equilibrum == true ) {
		for ( int i = 0; i < bandGapArray.size(); i++) {
			phip(i) = bandGapArray[i] - phin(i);
		}
	}
	// std::cout << materialList.size() << std::endl;
    int k=0;
    for ( int i=0; i < materialList.size(); i++) {
    	// std::cout << nyList[i] << std::endl;
    	// std::cout << typeArray[i] << std::endl;
        for ( int j=0; j < nyList[i]; j++) {
            if (typeList[i] == Semiconductor) {
                density(k) = materialList[i].chargeDensity((double)phin(k), (double)phip(k)); //TODO: why here used to be k+1, then "carrierDensity[0]=carrierDensity[1]; nyList has +1 on the bottom"
//                std::cout << (double)phin(k) << ", " << (double)phip(k) << std::endl;
//                std::cout << (double)density(k) << std::endl;
            } else if (typeList[i]==Dielectric) {
             };
            k++;
        };
    };

    return ( density );
}

sp_mat Device1D::qCMatFunct(mat phin, mat phip, bool Equilibrum) {
	if ( phin.n_elem != sumPoint || ( phip.n_elem != sumPoint ) )
		std::cerr << "Error: input wrong phin/phip size!" << std::endl;
    sp_mat cp(sumPoint, sumPoint); // bug 11/16: speye -> zeros

	if (Equilibrum == true ) {
		for ( int i = 0; i < bandGapArray.size(); i++) {
			phip(i) = bandGapArray[i] - phin(i);
		}
	}

    int k=0;
    for ( int i=0; i < materialList.size(); i++) {
        for ( int j=0; j < nyList[i]; j++) {
            if (typeList[i] == Semiconductor) {
                cp(k + 1, k + 1) = materialList[i].quantumCapa((double)phin(k), (double)phip(k)); //TODO: why here used to be k+1, then "cp(0, 0) = cp(1, 1);"
            } else if (typeList[i]==Dielectric) {
             };
            k++;
        };
    };
//     std::cout << cp << std::endl;
    return ( cp );
}


