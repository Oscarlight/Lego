/*
 * Poisson2D.h
 *
 *  Created on: Nov 4, 2015
 *      Author: oscar
 */

#ifndef SRC_POISSON2D_H_
#define SRC_POISSON2D_H_


#include "Device2D.h"
#include <map>
#include "Poisson1D.h" // TODO: should I include the parent class
#include "InOut2D.h"
/**
 * The rule of constructor call in derived class:
    1. Memory for cDerived is set aside (enough for both the Base and Derived portions).
    2. The appropriate Derived constructor is called
    3. The Base object is constructed first using the appropriate Base constructor !!!
    4. The initialization list initializes variables
    5. The body of the constructor executes
    6. Control is returned to the caller

    Inherited from Poisson1D since using setFLnArray, setFLpArray, bCArrayFunct,
      rangeByNum and rangeByStep from Poisson1D

    add default constructor in Poisson1D (base class constructor)
 */

class Poisson2D : public Poisson1D {
private:
	Device2D dev2D;
	// key: the block index; value Vtg, Wft, Vbg, Wfb
	std::map<int, std::array<double, 4>> biasMap;
	mat fLnArray; // like -Vds array in old version
	mat fLpArray;
	mat bCArray; // boundaryConditionArray
	//
	mat cDA; // chargeDensityArray
	mat phin;
	mat phip;
	//
	mat condBand;
	mat valeBand;

public:
	Poisson2D(Device2D _dev2D);
	virtual ~Poisson2D();
	/* Map(<index of the 1D block>, <fLnArray for a 1D block>) */
	void setFLnArray_2D(std::map<int, std::vector<double>> fLnMap); // if in equilibrium, only need Efn = -V
	void setFLpArray_2D(std::map<int, std::vector<double>> fLpMap);
	void setGateBias_2D(std::map<int, std::array<double, 4>> _biasMap);
	void runPoisson2D(double vTolerance, double _chargeTolerance, double magicNumber, bool Equilibrum);

	mat getBCArray_2D();
	mat getFLnArray_2D();
	mat getFLpArray_2D();
	mat getCondBand_2D();
	mat getValeBand_2D();
	mat getPhin_2D();
	mat getPhip_2D();
	mat getCDA_2D();

	// Helpers
	std::map<int, std::vector<double>> createFermiLevelMap(std::vector<bool> connect2Drain, double Vds, std::vector<bool> connect2Source, double Vss);
	// assume each 1D block has the same y-direction meshing
	// TODO: it is not generalized yet, only good for pin-TFET or Thin-TFET
	std::vector<int> sourceDrain2DLayerIndex(std::vector<bool> connect2Drain, std::vector<bool> connect2Source);

protected:
	mat bCArrayFunct_2D(Device2D dev2D, std::map<int, std::array<double, 4>> biasMap);
	mat setPhip_2D(mat phin, bool Equilibrum);
};

#endif /* SRC_POISSON2D_H_ */
