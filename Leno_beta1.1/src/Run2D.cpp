/*
 * Run2D.cpp
 *
 *  Created on: Nov 5, 2015
 *      Author: oscar
 */

#include "Run2D.h"

Run2D::Run2D() {
	// TODO Auto-generated constructor stub

}

Run2D::~Run2D() {
	// TODO Auto-generated destructor stub
}

std::map<std::string, Material> Run2D::readInput(int argc, char** argv) {
	GetPot commandLine(argc, argv);
	std::string Device2DFile;
	std::string MaterialFile;
    if( commandLine.search(2, "-2d", "--file") && commandLine.search(2, "-m", "--file") ) {
    	Device2DFile = commandLine.follow("", 2 , "-2d", "--file");
    	MaterialFile = commandLine.follow("", 2, "-m", "--file");
    }
    else {
    	cerr << "Error: input file missing. Use \" Leno -2d device2DFile -m materialFile \" " << endl ;
        exit(0);
    }

    // Create io
    io2D.readDevice2D(Device2DFile);
    io2D.readMaterial(MaterialFile);

	// create Materials, with default T = 300 K
	Params param;
	std::map<std::string, Material> matLib;
	int i = 0;
	for (std::map<std::string, std::vector<double>>::iterator it = io2D.matMap.begin();
			it  != io2D.matMap.end(); ++it) {
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

	std::cout << " ** material library is completed." << std::endl;
	return ( matLib );
}

Device2D Run2D::createDevice2D(int argc, char** argv) {
	std::map<std::string, Material> matLib = readInput(argc, argv);
	std::vector<Device1D> dev1DList;
	Device2D dev2D;
	int accu = 0; // accumalated layer num

	// construct all device1D
	for (int i = 0; i < io2D.blockNum; i++) { // for each block
		dev1DList.push_back(createDevice1D(io2D, matLib, accu, i));
		accu = accu + io2D.layerNumList[i]; // update accu
	}

	std::cout << " **** dev1DList is completed with size = " << dev1DList.size() << std::endl;

	// construct device2D
	dev2D.startWith_2D(dev1DList[0], io2D.blockLength[0], io2D.blockPoint[0], io2D.leftBT);

	for (int i = 1; i < io2D.blockNum -1; i++) {
		dev2D.add_2D(dev1DList[i], io2D.blockLength[i], io2D.blockPoint[i]);
	}
	dev2D.endWith_2D(dev1DList[io2D.blockNum - 1], io2D.blockLength[io2D.blockNum - 1],
			io2D.blockPoint[io2D.blockNum - 1], io2D.rightBT);
	dev2D.matrixDiff_2D();

	std::cout << " **** Device2D is completed." << std::endl;
	return ( dev2D );
}

std::vector<ExtractData> Run2D::getData2D() {
	return (data2DPerBias);
}

std::vector< std::vector< std::vector<double> > > Run2D::getBand2D() {
	return (band2DPerBias);
}

InOut2D Run2D::getIO2D() {
	return (io2D);
}

Device1D Run2D::createDevice1D(InOut2D io2D, std::map<std::string, Material> matLib, int accu, int i){
	Device1D dev1D;
    // construct device1D
    dev1D.startWith(matLib.at(io2D.layerName[accu]), io2D.layerThickness[accu], io2D.layerPoint[accu], io2D.botBTList[i]);
//    std::cout << io2D.layerName[accu] << "," << io2D.layerThickness[accu] << "," << io2D.layerPoint[accu]
//			  << "," << io2D.botBTList[i] << std::endl;
	for (int j = 1; j < io2D.layerNumList[i] - 1; j++) { // for each layer in the block
		dev1D.add(matLib.at(io2D.layerName[accu + j]), io2D.layerThickness[accu + j], io2D.layerPoint[accu + j]);
//		std::cout << io2D.layerName[accu + j] << "," << io2D.layerThickness[accu + j] << "," << io2D.layerPoint[accu + j] << std::endl;
	}
	dev1D.endWith(matLib.at(io2D.layerName[accu + io2D.layerNumList[i] -1]), io2D.layerThickness[accu + io2D.layerNumList[i] -1],
			io2D.layerPoint[accu + io2D.layerNumList[i] -1], io2D.topBTList[i]);

	return ( dev1D );
}

void Run2D::runPoisson2D(int argc, char** argv) {
	Device2D dev2D = createDevice2D(argc, argv);
	Poisson2D p2D(dev2D);
    std::vector<double> vtgArray =p2D.rangeByStep(io2D.vtgArray[0], io2D.vtgArray[1], io2D.vtgArray[2]);
    std::vector<double> vbgArray =p2D.rangeByStep(io2D.vbgArray[0], io2D.vbgArray[1], io2D.vbgArray[2]);
    std::vector<double> vdsArray =p2D.rangeByStep(io2D.vdsArray[0], io2D.vdsArray[1], io2D.vdsArray[2]);
    std::vector<double> vssArray =p2D.rangeByStep(io2D.vssArray[0], io2D.vssArray[1], io2D.vssArray[2]);

    std::vector<int> sdIndex = p2D.sourceDrain2DLayerIndex(io2D.connect2Drain, io2D.connect2Source);
    std::cout << "drain index: " << sdIndex[0] << " , source index: " << sdIndex[1] << std::endl;

    std::cout << " ***** Poisson2D is started." << std::endl;
    Capacitance capa;
    for (double Vbg : vbgArray) {
    	for (double Vss : vssArray) {
    		for (double Vds : vdsArray) {
    			for (double Vtg : vtgArray) {

    				if (io2D.terminalConnect[0] == true)  // Vbg = Vtg always
    					Vbg = Vtg;

    				std::cout << "Vbg = " << Vbg << " V; " << "Vtg = " << Vtg << " V; "
    						<< "Vds = " << Vds << " V; " << "Vss = " << Vss << std::endl;

    				p2D.setGateBias_2D(io2D.gateBiasMap(Vtg, Vbg));

//        			std::vector<double> fLnZero(dev2D.getUnitSize(), 0);   // 0 for one slice
//        			std::vector<double> fLn(dev2D.getUnitSize(), 0);   // 0 for one slice
//        			fLn.at(1 + io2D.layerPoint[io2D.layerNumList[1]] + 1 + io2D.layerPoint[io2D.layerNumList[1]+2]) = - Vds;
//        			std::map<int, std::vector<double>> fLnMap;
//        			fLnMap.insert(std::pair<int, std::vector<double>>(0, fLnZero));
//        			fLnMap.insert(std::pair<int, std::vector<double>>(1, fLn));
//        			fLnMap.insert(std::pair<int, std::vector<double>>(2, fLn));
//        			fLnMap.insert(std::pair<int, std::vector<double>>(3, fLn));
//        			p2D.setFLnArray_2D(fLnMap);
//        			p2D.setFLpArray_2D(fLnMap); //

    				/* -Vds: Caution: the negative sign */
    				p2D.setFLnArray_2D(p2D.createFermiLevelMap(io2D.connect2Drain, -Vds, io2D.connect2Source, Vss));
    				p2D.setFLpArray_2D(p2D.createFermiLevelMap(io2D.connect2Drain, -Vds, io2D.connect2Source, Vss)); // when equilibrum, Flp = Fln

    				p2D.runPoisson2D(io2D.voltageErr, io2D.carrierConcenErr, io2D.magicNum,  true);
    				std::cout << "2D Poisson Finished." << std::endl;

    				// Extract Data
    				ExtractData data;
    				data.bandAndCharge2D(p2D, dev2D); // first read into band info

    				// Store the bandAlignment
    				std::vector< std::vector<double> > tempBand;
    				tempBand.push_back(data.mapRowSlide(data.cBMap, sdIndex[0])); // top for ThinTFET
    				tempBand.push_back(data.mapRowSlide(data.vBMap, sdIndex[0]));
    				if (sdIndex[0] != sdIndex[1]) { // source and drain are not at the same layer
    					tempBand.push_back(data.mapRowSlide(data.cBMap, sdIndex[1])); // bottom for ThinTFET
    					tempBand.push_back(data.mapRowSlide(data.vBMap, sdIndex[1]));
    				}

    				band2DPerBias.push_back(tempBand);


    				// read top gate charge into capa for capacitance calculation
    				capa.readCharge(Vtg, Vbg, Vds, Vss, data.topGateCharge2D(dev2D, io2D.topGateArea)); // read in all the data needed for capa

    			}
    		}
    	}
    	if (io2D.terminalConnect[0] == true)  // Vbg = Vtg always
    		break;
    }

    capa.calCggCgsCgd(vtgArray, vbgArray, vdsArray, vssArray, io2D.terminalConnect); // calculate capacitance

    io2D.writeCapaMap(io2D.devName+"_"+io2D.userCom+"_Cgs", capa.getCgsMap());
    io2D.writeCapaMap(io2D.devName+"_"+io2D.userCom+"_Cgd", capa.getCgdMap());
    io2D.writeCapaMap(io2D.devName+"_"+io2D.userCom+"_Cgg", capa.getCggMap());
}
