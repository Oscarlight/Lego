/*
 * run1D.cpp
 *
 *  Created on: Nov 4, 2015
 *      Author: oscar
 */

#include "Run1D.h"

Run1D::Run1D() {
}

Run1D::~Run1D() {
}

std::map<std::string, Material> Run1D::readInput(int argc, char** argv) {
	GetPot  commandLine(argc, argv);
	std::string DeviceFile;
	std::string MaterialFile;

	if( commandLine.search(2, "-d", "--file") && commandLine.search(2, "-m", "--file") ) {
		DeviceFile = commandLine.follow("", 2 , "-d", "--file");
		MaterialFile = commandLine.follow("", 2, "-m", "--file");
	}
	else {
		cerr << "Error: input file missing. Use \" Leno -d deviceFile -m materialFile \" " << endl ;
		exit(0);
	}

	// Create io
	io.readDevice(DeviceFile);
	io.readMaterial(MaterialFile);
	// io.printInfo();

	// create Materials, with default T = 300 K
	Params param;
	std::map<std::string, Material> matLib;
	int i = 0;
	for (std::map<std::string, std::vector<double>>::iterator it = io.matMap.begin();
			it  != io.matMap.end(); ++it) {
		Material t(i, (int)round(it->second[0]), it->first.c_str(), it->second[1], it->second[2],
				it->second[3], it->second[4], (int)round(it->second[5]), it->second[6],
				it->second[7], it->second[8], it->second[9], it->second[10],
				it->second[11], it->second[12]);
		//		std::cout << i << "," << (int)round(it->second[0]) << "," << it->first.c_str() << "," << it->second[1] << ","
		//				  << it->second[2] << "," << it->second[3] << "," << it->second[4] << "," << (int)round(it->second[5]) << "," << it->second[6]
		//				  << "," << it->second[7] << "," << it->second[8] << "," << it->second[9] << "," << it->second[10] << "," << it->second[11] << ","
		// 		          << it->second[12] << std::endl;
		matLib.insert(std::pair<std::string, Material>(it->first, t));
	}

	return ( matLib );
}

Device1D Run1D::createDevice1D(int argc, char** argv) {
	Device1D dev1D;
	std::map<std::string, Material> matLib = readInput(argc, argv);
	// construct device1D
	dev1D.startWith(matLib.at(io.layerName[0]), io.layerThickness[0], io.layerPoint[0], io.botBT);
	// std::cout << io.layerName[0] << "," << io.layerThickness[0] << "," << io.layerPoint[0] << std::endl;
	for (int i = 1; i < io.layerNum-1; i++) {
		dev1D.add(matLib.at(io.layerName[i]), io.layerThickness[i], io.layerPoint[i]);
		// std::cout << io.layerName[i] << "," << io.layerThickness[i] << "," << io.layerPoint[i] << std::endl;
	}
	dev1D.endWith(matLib.at(io.layerName[io.layerNum-1]), io.layerThickness[io.layerNum-1],
			io.layerPoint[io.layerNum-1], io.topBT);
	// std::cout << io.layerName[io.layerNum-1] << "," << io.layerThickness[io.layerNum-1] << "," << io.layerPoint[io.layerNum-1] << std::endl;
	dev1D.matrixDiff();
	return ( dev1D );
}


void Run1D::runPoisson1D(int argc, char** argv){
	// calculate Poisson1D
	Device1D dev1D = createDevice1D(argc, argv);
	Poisson1D p1D(dev1D);
	std::vector<double> vtgArray =p1D.rangeByStep(io.vtgArray[0], io.vtgArray[1], io.vtgArray[2]);
	std::vector<double> vbgArray =p1D.rangeByStep(io.vbgArray[0], io.vbgArray[1], io.vbgArray[2]);
	std::vector<double> vdsArray =p1D.rangeByStep(io.vdsArray[0], io.vdsArray[1], io.vdsArray[2]);
	std::vector<double> vssArray =p1D.rangeByStep(io.vssArray[0], io.vssArray[1], io.vssArray[2]);

	for (double Vbg : vbgArray) {
		for (double Vss : vssArray) {
			for (double Vds : vdsArray) {
				Tunnelling t;
				for (double Vtg : vtgArray) {
					std::cout << "Vbg = " << Vbg << " V; " << "Vss = " << Vss << " V; "
							<< "Vtg = " << Vtg << " V; " << "Vds = " << Vds << " V;" << std::endl;

					p1D.setGateBias(Vtg, io.topWF, Vbg, io.botWF);
					/* Old manual mode to set Fermi level: */
//					std::vector<double> fLn1(dev1D.getSumPoint(), 0);
//					fLn1.at(1 + io.layerPoint[0] + 1 + io.layerPoint[2]) = - Vds;

					/* Updated method (02/18/2016) */
					std::vector<double> fLn = p1D.createFLnArray(io.connect2Drain, Vds,
							io.connect2Source, Vss);

					/*
					for (int i = 0; i < fLn1.size(); i++) {
						if (fLn1[i] != fLn[i]) {
							std::cout << "At " << i << " fLn1 = "
									<< fLn1[i] << " fLn = " <<
									fLn[i] << std::endl;
						}
					}
					*/

					p1D.setFLnArray(fLn);
					p1D.setFLpArray(fLn); // in Equilibrum Fln = Flp

					p1D.runPoisson1D(io.voltageErr, io.carrierConcenErr, io.magicNum,  true);
					std::cout << "Poisson Finished." << std::endl;

					// extract data
					ExtractData data;
					data.bandAndCharge(p1D, dev1D);
					if (io.calcuIV) {
						// current
						data.bASemiOnly();
						double iJ = t.interTunnel2DOFP(data.cBSemi[1],
								data.vBSemi[1], data.cBSemi[0], data.vBSemi[0],
								data.fLnSemi[0] - data.fLnSemi[1], 1E-3);
						double lJ = t.likeTunnel2DOFP(data.cBSemi[1],
								data.vBSemi[1], data.cBSemi[0], data.vBSemi[0],
								data.fLnSemi[0] - data.fLnSemi[1], 1E-3);

						// Store current in InputOut
						std::vector<double> current;
						current.push_back(iJ);
						current.push_back(iJ + lJ);
						io.storeIV(Vtg, Vbg, Vds, current);
					}

					// Store band data
					std::vector<std::vector<double>> tempBand;
					tempBand.push_back(data.x);
					tempBand.push_back(data.cB);
					tempBand.push_back(data.vB);
					tempBand.push_back(data.fLn);
					tempBand.push_back(data.fLp);
					tempBand.push_back(data.chargeDensity);
					tempBand.push_back(data.mobileElectronDensity);
					bandPerBias.push_back(tempBand);

				}
			}
		}
	}

	if (io.saveBA == true) {
		io.writeBAandCharge(io.devName+"_"+io.userCom+"_Bands_and_Charges", vbgArray, vdsArray, vtgArray, bandPerBias);
	}

}

InOut Run1D::getIO() {
	return ( io );
}

std::vector<std::vector<std::vector<double>>> Run1D::getBand(){
	return ( bandPerBias );
}
