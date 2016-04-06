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

#include "ExtractData.h"
#include "Poisson1D.h"
#include "Tunnelling.h"
#include "InOut.h"
#include "InOut2D.h"
#include "Run1D.h"
#include "Run2D.h"
/**
 * matplotlibcpp has to be include on the top level since it will cause conflicts/
 * multiple definision error when include in ExtraData class
 */
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
int main( int argc, char** argv ) {

	GetPot commandLine(argc, argv);
	if (commandLine.search(1, "-d")) {
		std::cout << "Enter 1D Mode" << std::endl;
		Run1D r1d;
		r1d.runPoisson1D(argc, argv);

		/* plotting and saving */
		InOut io = r1d.getIO();
		std::vector<std::vector<std::vector<double>>> dataArray = r1d.getBand();
		for (int i = 0; i < dataArray.size(); i++) {

			//band alignment plot
			if (io.plotBA == true) {
				plt::plot(dataArray[i][0], dataArray[i][1]);
				plt::plot(dataArray[i][0], dataArray[i][2]);
				plt::plot(dataArray[i][0], dataArray[i][3], "r--");
				plt::plot(dataArray[i][0], dataArray[i][4], "b--");
				plt::show();
			}
			//		if (io.savePl == true) {
			//			//TODO save throw runtime_err
			//			// plt::save("./plots/BA/BA.png");
			//		}
			//		if (io.plotIV == true) {
			//			//TODO Sometime when program is terminated in the middle, plots can not be closed.
			//			plt::plot(vtgArray, iArray);
			//			plt::show();
			//		}
		}
		//	if (io.savePl == true) {
		//		//TODO save throw runtime_err
		//		// plt::save("./plots/IV/I-Vds.png");
		//	}
		if (io.calcuIV == true) {
			io.writeIV((io.devName+"_"+io.userCom), io.getIVMap());
			std::cout << "--------End---------"<< std::endl;
		}

	} else if (commandLine.search(1, "-2d")) {
		std::cout << "Enter 2D Mode" << std::endl;
		Run2D r2d;
		r2d.runPoisson2D(argc, argv);

//		/* ploting and saving */
//		std::vector< std::vector< std::vector<double> > > bands = r2d.getBand2D();

//		for (int i = 0; i < bands.size(); i++) { // each bias
			//TODO: add a input to choose whether to plot this
//			for (int j = 0; j < bands[i].size(); j++) { // each bands
//			plt::plot(bands[i][j]);
//			}
//			plt::show();
//		}
//		plt::plot(tunnelWindows);
//		plt::show();
	}


}



