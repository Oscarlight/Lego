/*
 * InOut2D.cpp
 *
 *  Created on: Nov 12, 2015
 *      Author: oscar
 */

#include "InOut2D.h"

InOut2D::InOut2D() {
	// TODO Auto-generated constructor stub

}

InOut2D::~InOut2D() {
	// TODO Auto-generated destructor stub
}

void InOut2D::readDevice2D(std::string filename) {
	// open file, specify the comments are within /* */
	GetPot ifile(filename.c_str(), "/*", "*/");
	std::cout << "Device2D File is " << filename.c_str() << std::endl;
	//Device: general
	ifile.set_prefix("Device2D/"); // enter [Device2D] section
	devName = ifile("Device_Name","No name specified");
	std::cout << "Device Name is " << devName << std::endl;
	userCom = ifile("User_Comment", "No Comment");
	std::cout << "Comments: " << userCom << std::endl;
	blockNum = ifile("Block_Number", 0);
	std::cout << "Block Num is " << blockNum << std::endl;

	// Block length and points
	for (int i = 1; i <= blockNum; i++) {
		std::string str = "Block" + std::to_string(i) + "/";
		ifile.set_prefix(str.c_str());

		// layer number: how many layers within a block
		layerNumList.push_back(ifile("Layer_Number", 0));

		// set block length/points
		blockLength.push_back(ifile("Length", 0,0));
		blockPoint.push_back(ifile("Nx", 0));

		// boundary conditions
		if (ifile("Top_Boundary_Type", "") == "Dirichlet" ) {
			topBTList.push_back(Dirichlet);
			topWFList.push_back(ifile("Top_Gate_Workfunction", 0.0));
			topGateArea.push_back(i-1); // i - 1 because the first element is 0 internally
		} else if (ifile("Top_Boundary_Type", "") == "Neumann") {
			topBTList.push_back(Neumann);
			topWFList.push_back(0); // WF = 0 if Neumann
		} else
			std::cerr << "Error: top boundary type undefined." << std::endl;

		if (ifile("Bottom_Boundary_Type", "") == "Dirichlet" ) {
			botBTList.push_back(Dirichlet);
			botWFList.push_back(ifile("Bottom_Gate_Workfunction", 0.0));
		} else if (ifile("Bottom_Boundary_Type", "") == "Neumann"){
			botBTList.push_back(Neumann);
			botWFList.push_back(0); // WF = 0 if Neumann
		} else
			std::cerr << "Error: bottom boundary type undefined." << std::endl;
	}

	// Layers
	// for each block, look into each layers (two for-loop)
	for (int i = 1; i <= blockNum; i++) {
		for (int j = 1; j <= layerNumList[i-1]; j++ ) {
			std::string str = "Layer"+std::to_string(i) + std::to_string(j) +"/";
			ifile.set_prefix(str.c_str());
			layerName.push_back(ifile("Material", ""));
			layerThickness.push_back(ifile("Thickness",0.0));
			layerPoint.push_back(ifile("Ny",0));

			// Remember which layers are connect to drain/source
			(ifile("ConnectTo", "") == "Drain") ? connect2Drain.push_back(true) : connect2Drain.push_back(false);
			(ifile("ConnectTo", "") == "Source") ? connect2Source.push_back(true) : connect2Source.push_back(false);

			// first block need to define left boundary condition for each layer
			if (i == 1) {
				// left boundary condition
				if (ifile("Left_Boundary_Type", "") == "Dirichlet" ) {
					leftBT.push_back(Dirichlet);
				} else if (ifile("Left_Boundary_Type", "") == "Neumann"){
					leftBT.push_back(Neumann);
				} else
					std::cerr << "Error: left boundary type undefined." << std::endl;
			}

			// last block need to define right boundary condition for each layer
			if (i == blockNum) {
				// right boundary condition
				if (ifile("Right_Boundary_Type", "") == "Dirichlet" ) {
					rightBT.push_back(Dirichlet);
				} else if (ifile("Right_Boundary_Type", "") == "Neumann"){
					rightBT.push_back(Neumann);
				} else
					std::cerr << "Error: Right boundary type undefined." << std::endl;
			}
		}
	}

	// Read Bias
	ifile.set_prefix("Vtg/");
	vtgArray[0] = ifile("Vtg_From", 0.0);
	vtgArray[1] = ifile("Vtg_To", 0.0);
	vtgArray[2] = ifile("Vtg_Step", 0.0);
	ifile.set_prefix("Vbg/");
	vbgArray[0] = ifile("Vbg_From", 0.0);
	vbgArray[1] = ifile("Vbg_To", 0.0);
	vbgArray[2] = ifile("Vbg_Step", 0.0);
	ifile.set_prefix("Vds/");
	vdsArray[0] = ifile("Vds_From", 0.0);
	vdsArray[1] = ifile("Vds_To", 0.0);
	vdsArray[2] = ifile("Vds_Step", 0.0);
	ifile.set_prefix("Vss/");
	vssArray[0] = ifile("Vss_From", 0.0);
	vssArray[1] = ifile("Vss_To", 0.0);
	vssArray[2] = ifile("Vss_Step", 0.0);

	// Read Terminal Connection
	ifile.set_prefix("ConnectTerminal/");
	if (ifile("TopGate_and_BackGate", "") == "Yes")
		terminalConnect[0] = true;
	else if (ifile("TopGate_and_BackGate", "") == "No")
		terminalConnect[0] = false;
	else
		std::cerr << "Error: TopGate_and_BackGate is not defined." << std::endl;
	if (ifile("TopGate_and_Drain", "") == "Yes")
		terminalConnect[1] = true;
	else if (ifile("TopGate_and_Drain", "") == "No")
		terminalConnect[1] = false;
	else
		std::cerr << "Error: TopGate_and_Drain is not defined." << std::endl;
	if (ifile("BackGate_and_Drain", "") == "Yes")
		terminalConnect[2] = true;
	else if (ifile("BackGate_and_Drain", "") == "No")
		terminalConnect[2] = false;
	else
		std::cerr << "Error: BackGate_and_Drain is not defined." << std::endl;


	// Read Poisson Convergence Condition
	ifile.set_prefix("PoissonConvergence/");
	voltageErr = ifile("Voltage_Error", 0.0);
	carrierConcenErr = ifile("Carrier_Density_Error", 0.0) * 1E10; // in cm^-3
	magicNum = ifile("Magic_Number", 0.0);
}

void InOut2D::readMaterial(std::string filename) {
	// open file, specify the comments are within /* */
	GetPot ifile(filename.c_str(), "/*", "*/");
	std::cout << "Material File is " << filename.c_str() << std::endl;
	// get the materials device structure required
	for (int i = 0; i < layerName.size(); i++) {
		std::vector<double> v;
		std::string str = layerName[i] + "/";
		std::cout << "Getting " << layerName[i] << std::endl;
		ifile.set_prefix(str.c_str());

		if (ifile("Type","") == "Dielectric")
			v.push_back(Dielectric);
		else if (ifile("Type","") == "Semiconductor")
			v.push_back(Semiconductor);
		else
			std::cerr << "Error: Wrong type" << std::endl;

		v.push_back(ifile("er_y",0.0));
		v.push_back(ifile("er_x",0.0));
		v.push_back(ifile("Electron_Affinity",0.0));
		v.push_back(ifile("Bandgap",0.0));
		v.push_back(ifile("Dimension",0.0));
		v.push_back(ifile("Monolayer_Thickness",0.0));
		v.push_back(ifile("Nd",0.0));
		v.push_back(ifile("Na",0.0));
		// Feb 6th, 2017. Uncomment Nta and Ntd lines. (By Frank)
		v.push_back(ifile("Nta",0.0));
		v.push_back(ifile("Ntd",0.0));
		v.push_back(ifile("mc_eff",0.0));
		v.push_back(ifile("mv_eff",0.0));
		v.push_back(ifile("E_CNL",0.0));
		v.push_back(ifile("gv_of_Ec",0.0));
		v.push_back(ifile("gv_of_Ev",0.0));

		matMap.insert(std::pair<std::string, std::vector<double>>(layerName[i], v));
		// Map will eliminate duplication automatically

	}
	for (std::map<std::string, std::vector<double>>::iterator it = matMap.begin();
			it  != matMap.end(); ++it) {
		std::cout << "Material: " << it->first.c_str() << " is added." << std::endl;
	}

}

std::map<int, std::array<double, 4>> InOut2D::gateBiasMap(double Vtg, double Vbg) {
	// Construct BiasMap from bias arrays and Workfunction List
	// for the Block boundary condition = Neumann, Voltage = 0; Workfunction = 0;
	std::map<int, std::array<double, 4>> biasMap;
	for (int i = 0; i < blockNum; i++) {
		// add bias and workfunction
		std::array<double, 4> tempGateBias = {Vtg, topWFList[i], Vbg, botWFList[i]};
		biasMap.insert(std::pair<int, std::array<double, 4>>(i, tempGateBias));
	}
	return ( biasMap );
}

void InOut2D::writeCapaMap(std::string fileName, std::map<std::vector<double>, double> map) {
	std::ofstream myfile;
	myfile.open((fileName+".csv").c_str());
	myfile << "Vbg" << ", " << "Vtg for Cgd/Cgs (Vds for Cgg)" << ", " << "Vds for Cgd/Cgs (Vtg for Cgg)" << ", " << "Capacitance (fF/um)" << "\n";
	for (std::map<std::vector<double>, double>::iterator it = map.begin(); it != map.end(); ++it) {
		for (int i = 0; i < it->first.size(); i++) {
			// double to string method 2
			myfile << std::to_string(it->first[i]) << ", ";
		}
		myfile << std::to_string(it->second * 1E11)  << "\n"; // fF/um
	}
	myfile.close();
}

void InOut2D::writeGateEffMap(std::string fileName, std::map<std::vector<double>, std::vector<double>> map) {
	std::ofstream myfile;
	myfile.open((fileName+".csv").c_str());
	myfile << "Vbg" << ", " << "Vds" << ", " << "first tunnel window (eV)" << ", " << "last tunnel window (eV)" << ", "
	<< "efficiency around 0 (eV/V)" << ", " << "threshold voltage (V)" << "\n";
	for (std::map<std::vector<double>, std::vector<double>>::iterator it = map.begin(); it != map.end(); ++it) {
		for (int i = 0; i < it->first.size(); i++) {
			myfile << std::to_string(it->first[i]) << ", ";
		}
		for (int i = 0; i < it->second.size(); i++) {
			myfile << std::to_string(it->second[i])<< ", ";
		}
		myfile << "\n";
	}
	myfile.close();
}

/** write in a form of matrix over the span of the whole device,
 * normally only used for signal bias simulation */
void InOut2D::write2DData(std::string fileName,
		std::map<int, std::vector<std::vector<double> > > dataMap) {
	std::ofstream myfile;
	myfile.open((fileName+".csv").c_str());

	int unitSize = dataMap.at(0)[0].size();
	// std::cout << unitSize << std::endl;
	for (int k = 0; k < unitSize; k++) { // point
		for (int i = 0; i < dataMap.size(); i++) { // block
			for (int j = 0; j < dataMap.at(i).size(); j++) {// slice
				myfile << dataMap.at(i)[j][k] << ", ";
			}
		}
		myfile << "\n";
	}
	myfile.close();
	std::cout << fileName << " is exported." <<  std::endl;
}


void InOut2D::writeByLayerinSemi(std::string fileName, std::vector<double> vtgArray,
		std::vector<double> vdsArray, std::vector<double> vbgArray,
		std::vector< std::vector< std::vector<double> > > data2DPerBias) {

	std::ofstream myfile;
	myfile.open((fileName+".csv").c_str());
	myfile << "Vbg" << ", " << "Vds" << ", " << "Vtg" << ", " << "index" << ", " << "data"  << "\n";

	// in consistent with the for-loop when calculating 2D Poisson in Run2D.cpp
	int accu = 0;
	int numOfBands = data2DPerBias[0].size();
	int numOfPointInBand = data2DPerBias[0][0].size();
	for (int i = 0; i < vbgArray.size(); i++) {
		for (int j = 0; j < vdsArray.size(); j++) {
			for (int k = 0; k < vtgArray.size(); k++) {

				for (int m = 0; m < numOfBands; m++) {
					myfile << vbgArray[i] << ", " << vdsArray[j] << ", " << vtgArray[k] << ", " << m << ", ";
					for (int n = 0; n < numOfPointInBand; n++) {
						myfile << data2DPerBias[accu][m][n] << ", ";
					}
					myfile << "\n";
				}
				accu++;
			}
		}
	}

	myfile.close();
	std::cout << fileName << ".csv" << " is exported." <<  std::endl;
}


/* static */
void InOut2D::printMatToText(const char* filename, mat matrix)
{
    std::ofstream MyFile;
    MyFile.open(filename);
    MyFile << matrix;
    MyFile.close();
}
