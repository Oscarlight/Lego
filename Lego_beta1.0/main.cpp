/*
 * main.cpp
 *
 *  Created on: Oct 11, 2015
 *      Author: oscar
 */

/**
 * test Poisson 1D
 */
#include <stdio.h>
#include "Device1D.h"
#include "Poisson1D.h"
int main( int argc, const char* argv[] )
{
	// create Materials, with default T = 300 K
	double t2D = 0.6; // nm
	Params param;
	Material   WSe2(0, param.Semiconductor,  "WSe2", 1,     6,   4, 1.3, 2, t2D, 0, 0, 0.3, 0.4, 2, 2);
    Material  SnSe2(1, param.Semiconductor, "SnSe2", 3.9,   6, 5.1, 1,   2, t2D, 0, 0, 0.3, 0.4, 2, 2);
    Material   SiO2(2, param.Dielectric,     "SiO2", 3.9, 3.9,   1, 9,   3, 0,   0, 0,   0,   0, 0, 0);
    Material vdWGap(2, param.Dielectric,   "vdWGap", 1, 	  1, 4, 2,   3, 0,   0, 0,   0,   0, 0, 0);

    // create device1D
    Device1D dev1D;
    dev1D.startWith(SiO2, 1, 100, param.Dielectric);
    dev1D.add(WSe2, 0, 0);
    dev1D.add(vdWGap, 0.35, 35);
    dev1D.add(SnSe2, 0, 0);
    dev1D.endWith(SiO2, 1, 100, param.Dielectric);

    // calculate Poisson1D
    Poisson1D p1D(dev1D);
    p1D.setGateBias(0.2, 5.05, 0, 5.65);
    std::vector<double> fLn(dev1D.getSumPoint(), 0);
    fLn.at(137) = - 0.2; // + 0.2 V
    p1D.setFLnArray(fLn);
    p1D.runPoisson1D(0.001, 1E10, 1,  true);

}



