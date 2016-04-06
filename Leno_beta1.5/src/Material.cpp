/*
 * Material.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: oscar
 */

#include "Material.h"
Material::Material() {
	// std::cout << "Default constructor of Material is called." << std::endl;
}
/**
 * Constructor using default T
 */
Material::Material(int _index, int _type, const char *pname, double _dielectricConstant, double _dieleIP, double _electronAffinity, double _bandGap, int _dimension, \
				  double _t2D, double _Nd, double _Na, double _mcEff, double _mvEff, int _gvC, int _gvV) : \
				  index(_index), type(_type),name(pname), dieleY(_dielectricConstant), dieleX(_dieleIP), \
				  electronAffinity(_electronAffinity),bandGap(_bandGap),dimension(_dimension),
				  gvC(_gvC),gvV(_gvV) {
	    	std::cout << "Initialize Material " << pname  << std::endl;

	    	if (dimension == 2 && _t2D == 0)
	    		std::cerr << "Error: no native thickness defined for 2D material:" << pname  << std::endl;

	    	// Unit convert from natural unit to Lego unit
	    	t2D = _t2D * cmNL(1E-7); // origin: nm
	    	Nd = _Nd * 1 / ( cmNL(1E-7) * cmNL(1E-7) * cmNL(1E-7) ); // origin: 1/nm^3
	    	Na = _Na * 1 / ( cmNL(1E-7) * cmNL(1E-7) * cmNL(1E-7) );
	    	Nta = 0; // Acceptor-like trap density (same unit as Nd)
	    	Ntd = 0;
//	    	Nta = _Nta * 1 / ( cmNL(1E-7) * cmNL(1E-7) * cmNL(1E-7) ); // origin: 1/nm^3
//	    	Ntd = _Ntd * 1 / ( cmNL(1E-7) * cmNL(1E-7) * cmNL(1E-7) );
	        mcEff=_mcEff*m0;
	        mvEff=_mvEff*m0;
}

Material::Material(double _T, int _index, int _type, const char *pname, double _dielectricConstant, \
				  double _dieleIP, double _electronAffinity, double _bandGap, int _dimension, \
				  double _t2D, double _Nd, double _Na, double _mcEff, double _mvEff, int _gvC, int _gvV) : Params(_T) {\
            Material(_index, _type, pname, _dielectricConstant, _dieleIP, _electronAffinity, _bandGap, _dimension, _t2D, \
					         _Nd, _Na, _mcEff, _mvEff, _gvC, _gvV);
		    std::cout << "Initialize Material " << pname << " with T = " << _T << std::endl;
}


Material::~Material() {}

/**
 * return the electronDenstiy given phin (i.e. Ec - Efn)
 */
double Material::electronDensity(double phin) {
	double electronDensity=0;
	switch (dimension) {
	case 2:
	{
		electronDensity=gvC*mcEff*Vth*ChargeQ/(Pi*SQUARE(Planckba))*log1p(exp(-phin/Vth)) / t2D;
	}
	break;
	case 3:
	{
		//TODO: IMPORTANT!! ALWAYS get Nta = 0, do make clean, problem solved!!
		//TODO: add Nta options in the input file. tried but ExtractData will not run correctly
		electronDensity=2*gvC*pow( (mcEff*Vth*ChargeQ/(2*Pi*SQUARE(Planckba))), 1.5) *
				fermiIntegralHalf(-phin/Vth) + 1600 * exp(-(phin - 0.14)/(1 * Vth))*Vth*ChargeQ*log1p(exp((phin - 0.14)/(1 * Vth)));
//		electronDensity=2*gvC*pow( (mcEff*Vth*ChargeQ/(2*Pi*SQUARE(Planckba))), 1.5) * \
//						fermiIntegralHalf(-phin/Vth);
//		std::cout << " electronDensity = " << electronDensity << std::endl;
//		std::cout << " fermi Integral half = " << fermiIntegralHalf(-phin/Vth) << std::endl;
//		std::cout << " Nc = " << 2*gvC*pow( (mcEff*Vth*ChargeQ/(2*Pi*SQUARE(Planckba))), 1.5) << std::endl;
//		std::cout << " Nta_all = " << Nta*exp(-phin/Vth)*Vth*ChargeQ*log1p(exp(phin/Vth)) << std::endl;
//		std::cout << " Nta = " << Nta << std::endl;
//		std::cout << " exp(-phin/Vth) = " << exp(-phin/Vth) << std::endl;
//		std::cout << " Vth*ChargeQ = " << Vth*ChargeQ << std::endl;
//		std::cout << " log1p(exp(phin/Vth)) =  " << log1p(exp(phin/Vth)) << std::endl;
//		std::cout << " phin / Vtg = " << phin/Vth << std::endl;
	}
	break;
	default:
	{
		std::cerr << " Error: You can only choose either 2D or 3D for carrierDensity. " << std::endl;
	}
	break;
	}
	return ( electronDensity );
}

/**
 * return the electronDenstiy given phin (i.e. Ec - Efn)
 */
double Material::mobileElectronDensity(double phin) {
	double electronDensity=0;
	switch (dimension) {
	case 2:
	{
		electronDensity=gvC*mcEff*Vth*ChargeQ/(Pi*SQUARE(Planckba))*log1p(exp(-phin/Vth)) / t2D;
	}
	break;
	case 3:
	{
		//TODO: They are the eletron in the Conduction band, I didn't consider the subgap trap DOS first
		// so the orginal eletronDensity function has to include all the trap DOS, which are not going to
		// contribute to the mobile charges
		electronDensity=2*gvC*pow( (mcEff*Vth*ChargeQ/(2*Pi*SQUARE(Planckba))), 1.5) *
				fermiIntegralHalf(-phin/Vth);
	}
	break;
	default:
	{
		std::cerr << " Error: You can only choose either 2D or 3D for carrierDensity. " << std::endl;
	}
	break;
	}
	return ( electronDensity );
}


/**
 * return the holeDensity given phip (i.e. Efp - Ev)
 * 2D in 1/m2; 3D in 1/m3
 */
double Material::holeDensity(double phip) {
	double holeDensity=0;
	switch (dimension) {
	case 2:
	{
		holeDensity = gvV*mvEff*Vth*ChargeQ/(Pi*SQUARE(Planckba))*log1p(exp(-phip/Vth)) / t2D;
	}
	break;
	case 3:
	{
		//TODO: add Ntd options in the input file. tried but ExtractData will not run correctly
        holeDensity=2*gvV*pow( (mvEff*Vth*ChargeQ/(2*Pi*SQUARE(Planckba))), 1.5) *
        		fermiIntegralHalf(-phip/Vth) + Ntd*exp(-phip/Vth)*Vth*ChargeQ*log1p(exp(phip/Vth));
//        holeDensity=2*gvV*pow( (mvEff*Vth*ChargeQ/(2*Pi*SQUARE(Planckba))), 1.5) * \
//         		fermiIntegralHalf(-phip/Vth);
	}
	break;
	default:
	{
		std::cerr << " Error: You can only choose either 2D or 3D for carrierDensity. " << std::endl;
	}
	break;
	}
	return ( holeDensity );
}

/**
 * return the total charge Density for semiconductor
 */
double Material::chargeDensity(double phin, double phip)
{
    return ( Nd - Na + holeDensity(phip) - electronDensity(phin) );
};

/**
 * return the total (fixed) charge density for insulator
 */
double Material::fixChargeDensity()
{
	if (type != Dielectric) {
		std::cerr << " Error: fixChargeDensity() method is only defined for insulator" << std::endl;
	}
    return ( Nd - Na );
};


/**
 * return the quantum capacitance
 */
double Material::quantumCapa(double phin, double phip)
{
    double CqElectron = 0;
    double CqHole = 0;

    switch (dimension)
    {
        case 2:
        {
            CqElectron = - gvC*mcEff*ChargeQ/(Pi*SQUARE(Planckba))* 1 / ( 1 + exp(phin/Vth) ) * 1 / t2D;
            CqHole = gvV*mvEff*ChargeQ/(Pi*SQUARE(Planckba))* 1 / ( 1 + exp(phip/Vth) ) * 1 / t2D;
        }
            break;
        case 3:
        {
            CqElectron = 0; // do not consider Cq for 3 dimension for now
            CqHole = 0;
        }
            break;
        default:
        {
            // std::cout << " You can only choose either 2D or 3D for quantumCapa. " << std::endl;
        }
            break;
    }
    return ( CqHole - CqElectron );
};

