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
/**
 * matplotlibcpp has to be include on the top level since it will cause conflicts/
 * multiple definision error when include in ExtraData class
 */
#include "lib/matplotlib-cpp/matplotlibcpp.h"
namespace plt = matplotlibcpp;

int main( int argc, const char* argv[] )
{
	// create Materials, with default T = 300 K
	double t2D = 0.6; // nm
	Params param;
	Material   nGaN(0, param.Semiconductor,  "nGaN", 9,  9, 4.1, 3.4, 3, 0, 2e15, 0, 0.2, 1.5, 2, 4);
    Material   pGaN(1, param.Semiconductor,  "pGaN", 9,  9, 4.1, 3.4, 3, 0, 0, 7e16, 0.2, 1.5, 2, 4);

    // create device1D
    Device1D dev1D;
    dev1D.startWith(nGaN, 8E3, 100, param.Dirichlet);
    dev1D.endWith(pGaN, 400, 100, param.Dirichlet);
    dev1D.matrixDiff();

    // calculate Poisson1D
    Poisson1D p1D(dev1D);
    // p1D.setGateBias(0.2, 5.05, 0, 5.65);

    std::vector<double> fLn(dev1D.getSumPoint(), 2.7);
    for (int i = 0; i < 101; i++) {
    	fLn.at(i) = 0;
    }
    // fLn.at(137) = - 0.2; // + 0.2 V
    p1D.setFLnArray(fLn);
    p1D.setFLpArray(fLn); // in Equilibrum Fln = Flp
    p1D.runPoisson1D(0.001, 1E10, 1,  true);

    // extract data
    ExtractData data;
    data.bandAndCharge(p1D, dev1D);

    // plot
    plt::plot(data.x, data.cB);
    plt::plot(data.x, data.vB);
    plt::plot(data.x, data.fLn, "r--");
    plt::plot(data.x, data.fLp, "b--");
    plt::show();
    plt::plot(data.x, data.chargeDensity);
    plt::show();
}



