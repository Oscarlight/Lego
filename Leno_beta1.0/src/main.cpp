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
#include "ExtractData.h"
#include "Poisson1D.h"
#include "Transport.h"
/**
 * matplotlibcpp has to be include on the top level since it will cause conflicts/
 * multiple definision error when include in ExtraData class
 */
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

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
    dev1D.startWith(SiO2, 1, 100, param.Dirichlet);
    dev1D.add(SnSe2, 0, 0);
    //dev1D.add(WSe2, 0, 0);
    dev1D.add(vdWGap, 0.35, 35);
    //dev1D.add(SnSe2, 0, 0);
    dev1D.add(WSe2, 0, 0);
    dev1D.endWith(SiO2, 1, 100, param.Dirichlet);
    dev1D.matrixDiff();


    // calculate Poisson1D
    Poisson1D p1D(dev1D);
    std::vector<double> vtgArray =p1D.rangeByStep(-0.4, 0, 0.05);
    std::vector<double> vdsArray =p1D.rangeByStep(-0.4, -0.4, 0.05);
    Transport t;
    for (double Vtg : vtgArray) {
    	for (double Vds : vdsArray) {
    		// p1D.setGateBias(Vtg, 5.05, 0, 5.65);
    		p1D.setGateBias(Vtg, 5.35, 0, 4.8);
    		std::vector<double> fLn(dev1D.getSumPoint(), 0);
    		fLn.at(137) = - Vds;
    		p1D.setFLnArray(fLn);
    		p1D.setFLpArray(fLn); // in Equilibrum Fln = Flp
    		p1D.runPoisson1D(0.001, 1E10, 1,  true);

    		// extract data
    		ExtractData data;
    		data.bandAndCharge(p1D, dev1D);

    		//plot
    		    plt::plot(data.x, data.cB);
    		    plt::plot(data.x, data.vB);
    		    plt::plot(data.x, data.fLn, "r--");
    		    plt::plot(data.x, data.fLp, "b--");
    		    // plt::show();

    		// current
    		data.bASemiOnly();
    		t.interTunnel2DOFP(data.cBSemi[1], data.vBSemi[1], data.cBSemi[0], data.vBSemi[0], data.fLnSemi[0] - data.fLnSemi[1], 1E-3);
    	}
    }
    //plt::show();
    std::vector<double> iArray = t.getInterTunnelCurrent();
    plt::plot(vtgArray, iArray);
    plt::show();
}



