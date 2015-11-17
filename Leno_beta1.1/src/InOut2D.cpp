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
	for (int i = 1; i <= blockNum; i++) {
		for (int j = 1; j <= layerNumList[i-1]; j++ ) {
			std::string str = "Layer"+std::to_string(i) + std::to_string(j) +"/";
			ifile.set_prefix(str.c_str());
			layerName.push_back(ifile("Material", ""));
			layerThickness.push_back(ifile("Thickness",0.0));
			layerPoint.push_back(ifile("Ny",0));

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
		v.push_back(ifile("mc_eff",0.0));
		v.push_back(ifile("mv_eff",0.0));
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

void InOut2D::writeCB(std::string fileName,
		std::map<int, std::vector<std::vector<double> > > cBMap) {
	std::ofstream myfile;
	myfile.open((fileName+"_condBand.csv").c_str());

	int unitSize = cBMap.at(0)[0].size();
	// std::cout << unitSize << std::endl;
	for (int k = 0; k < unitSize; k++) { // point
		for (int i = 0; i < cBMap.size(); i++) { // block
			for (int j = 0; j < cBMap.at(i).size(); j++) {// slice
				// std::cout << cBMap.at(i).size() << std::endl; // # of slice in a block
				// std::cout <<  k << ", " << i << ", " << j << std::endl;
				myfile << cBMap.at(i)[j][k] << ", ";
			}
		}
		myfile << "\n";
	}
	myfile.close();
	std::cout << "Conduction band data: " << fileName << "_condBand.csv" << " is exported." <<  std::endl;
}

/* static */
void InOut2D::printMatToText(const char* filename, mat matrix)
{
    std::ofstream MyFile;
    MyFile.open(filename);
    MyFile << matrix;
    MyFile.close();
}
