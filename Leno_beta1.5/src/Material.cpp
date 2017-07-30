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
			double _t2D, double _Nd, double _Na, double _mcEff, double _mvEff, int _gvC, int _gvV, double _Nta, double _Ntd,  double _E_CNLn, double _E_CNLp, \
			double _Wta, double _Wtd, double _E_peakn, double _E_peakp) : \
				  index(_index), type(_type), name(pname), dieleY(_dielectricConstant), dieleX(_dieleIP), \
				  electronAffinity(_electronAffinity),bandGap(_bandGap),dimension(_dimension), E_CNLn(_E_CNLn), E_CNLp(_E_CNLp), \
				  gvC(_gvC),gvV(_gvV),Wta(_Wta),Wtd(_Wtd),E_peakn(_E_peakn),E_peakp(_E_peakp) {
	    	std::cout << "Initialize Material " << pname  << std::endl;

	    	if (dimension == 2 && _t2D == 0)
	    		std::cerr << "Error: no native thickness defined for 2D material:" << pname  << std::endl;

	    	// Unit convert from natural unit to Lego unit
	    	t2D = _t2D * cmNL(1E-7); // origin: nm
	    	Nd = _Nd * 1 / ( cmNL(1E-7) * cmNL(1E-7) * cmNL(1E-7) ); // origin: 1/nm^3
	    	Na = _Na * 1 / ( cmNL(1E-7) * cmNL(1E-7) * cmNL(1E-7) );
		// Feb 6th, 2017. Uncomment Nta and Ntd lines and comment Nta=0 and Ntd=0 lines. (By Frank)
	    	Nta = _Nta * 1 / ( cmNL(1E-7) * cmNL(1E-7) * cmNL(1E-7) ); // origin: 1/nm^3
	    	Ntd = _Ntd * 1 / ( cmNL(1E-7) * cmNL(1E-7) * cmNL(1E-7) );
	        mcEff=_mcEff*m0;
	        mvEff=_mvEff*m0;
		//std::cout<<"bandGap is "<< bandGap <<std::endl;
		//std::cout<<"Nta = "<<Nta<<std::endl;
		//std::cout<<"Ntd = "<<Ntd<<std::endl;
		//std::cout<<"mcEff = "<<mcEff<<std::endl;
		//std::cout<<"mvEff = "<<mvEff<<std::endl;
		//std::cout<<"gvC = "<<gvC<<std::endl;
		//std::cout<<"gvV = "<<gvV<<std::endl;
		//std::cout<<"E_CNLn = "<<E_CNLn<<std::endl;
		//std::cout<<"E_CNLp = "<<E_CNLp<<std::endl;
		//std::cout<<"Wta = "<<Wta<<std::endl;
		//std::cout<<"Wtd = "<<Wtd<<std::endl;
		//std::cout<<"E_peakn = "<<E_peakn<<std::endl;
		//std::cout<<"E_peakp = "<<E_peakp<<std::endl;
}

Material::Material(double _T, int _index, int _type, const char *pname, double _dielectricConstant, double _dieleIP, double _electronAffinity, double _bandGap, int _dimension, \
			double _t2D, double _Nd, double _Na, double _mcEff, double _mvEff, int _gvC, int _gvV, double _Nta, double _Ntd,  double _E_CNLn, double _E_CNLp, \
			double _Wta, double _Wtd, double _E_peakn, double _E_peakp) : Params(_T)  {
		Material(_index, _type, pname, _dielectricConstant, _dieleIP, _electronAffinity, _bandGap, _dimension, _t2D, \
					         _Nd, _Na, _mcEff, _mvEff, _gvC, _gvV, _Ntd, _Nta, _E_CNLn, _E_CNLp, _Wta, _Wtd, _E_peakn, _E_peakp);
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
		if (Nta-0.0 < 1e-14) {
			electronDensity=gvC*mcEff*Vth*ChargeQ/(Pi*SQUARE(Planckba))*log1p(exp(-phin/Vth)) / t2D;
		}
		// with Trap DOS (with E_CNL at the center of Eg),  1e12 1/(cm^2 * eV) -> 1 in Leno unit
		// Feb 14th, 2017. Change to the eq with E_CNL. (By Frank)
		else {
		// E_CNL here is the difference between 1/2 Eg and actuall E_CNL. If E_CNL is above 1/2 Eg, E_CNL is positive. Change 10 to Nta.
			electronDensity=gvC*mcEff*Vth*ChargeQ/(Pi*SQUARE(Planckba))*log1p(exp(-phin/Vth)) / t2D +
					Nta * Vth * log1p(exp((bandGap/2 - E_CNLn - phin)/Vth)) / t2D;
		}
		//std::cout<<"Nta = "<<Nta<<std::endl;
		//std::cout<<"Nta_all = "<<Nta*Vth*log1p(exp((bandGap/2-E_CNLn-phin)/Vth))/t2D <<std::endl;
	}
	break;
	case 3:
	{
		//TODO: IMPORTANT!! ALWAYS get Nta = 0, do make clean, problem solved!!
		//TODO: add Nta options in the input file. tried but ExtractData will not run correctly
		// for DRC2016 Nta = 1600
		// Feb 6th, 2017. Change 1600 to Nta. (By Frank)
		if (Nta-0.0 < 1e-14) {
			electronDensity=2*gvC*pow( (mcEff*Vth*ChargeQ/(2*Pi*SQUARE(Planckba))), 1.5) * fermiIntegralHalf(-phin/Vth);
		}
		else {
			// Modification: open port for input Wta (width of exponential distribution, default=0.026eV) and E_peakn (peak position of exponential distribution, default=0eV for electron)
			// IMPORTANT MODIFICATION: change the calculation of trap density. (Mar 8th, 2017)
			// The calculation is done by direct integral using MATLAB and fit the result into polynomial (5th order). The variety range is:
			// Wta = 0.026 to 1 eV.        E_peakn-phin = -1 to 1 eV.     E_peakn = Ec-Epeak
			double En = E_peakn - phin;   // En = Efn-Epeak
			if (fabs(En)>1)
				cout<<"WARNING: The difference between E_peakn and Efn is too large, which may cause inaccuracy."<<endl;
			double n_trap = - 0.0002727 + 0.08858*En + 0.9995*Wta + 0.004625*pow(En,2) + 2.973*En*Wta + 0.01115*pow(Wta,2) - 0.4243*pow(En,3) - 0.01151*pow(En,2)*Wta - 2.744*En*pow(Wta,2)
					- 0.01536*pow(Wta,3) - 0.007636*pow(En,4) - 2.475*pow(En,3)*Wta - 0.05922*pow(En,2)*pow(Wta,2) - 4.608*En*pow(Wta,3) + 0.06716*pow(Wta,4)
					+ 0.3877*pow(En,5) + 0.03381*pow(En,4)*Wta + 4.877*pow(En,3)*pow(Wta,2) + 0.05508*pow(En,2)*pow(Wta,3) + 5.089*En*pow(Wta,4) - 0.1897*pow(Wta,5);
			electronDensity = 2*gvC*pow( (mcEff*Vth*ChargeQ/(2*Pi*SQUARE(Planckba))), 1.5) * fermiIntegralHalf(-phin/Vth) + Nta * n_trap;
			//	+ Nta * exp(-(phin - E_peakn)/(nwn * Vth))*(nwn*Vth)*ChargeQ*log1p(exp((phin - E_peakn)/(nwn * Vth)));
		}
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
	//std::cout<<"electron Density = "<<electronDensity<<endl;
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
		if (Ntd-0.0 < 1e-14) {
			holeDensity = gvV*mvEff*Vth*ChargeQ/(Pi*SQUARE(Planckba))*log1p(exp(-phip/Vth)) / t2D;
		}
		// with Trap DOS (with E_CNL at the center of Eg),  1e12 1/(cm^2 * eV) -> 1 in Leno unit
		// Feb 14th, 2017. Change to the eq with E_CNL. (By Frank)
		else {
		// E_CNL here is the difference between 1/2 Eg and actuall E_CNL. If E_CNL is above 1/2 Eg, E_CNL is positive. Change 10 to Ntd.
			holeDensity = gvV*mvEff*Vth*ChargeQ/(Pi*SQUARE(Planckba))*log1p(exp(-phip/Vth)) / t2D +
					Ntd * Vth * log1p(exp((bandGap/2 + E_CNLp - phip)/Vth)) / t2D;
		}
		//std::cout<<"Ntd = "<<Ntd<<std::endl;
		//std::cout<<"Ntd_all = "<<Ntd*Vth*log1p(exp((bandGap/2+E_CNLp-phip)/Vth))/t2D << std::endl;
	}
	break;
	case 3:
	{
		//TODO: add Ntd options in the input file. tried but ExtractData will not run correctly
		// Feb 14th, 2017. Consider Ntd. (By Frank)
		if (Ntd-0.0 < 1e-14) {
		        holeDensity = 2*gvV*pow( (mvEff*Vth*ChargeQ/(2*Pi*SQUARE(Planckba))), 1.5) * fermiIntegralHalf(-phip/Vth);
		}
		else {
			// Modification: open port for input Wtd (width of exponential distribution, default=0.026) and E_peakp (peak position of exponential distribution, default=0 for hole)
			// IMPORTANT MODIFICATION: change the calculation of trap density. (Mar 8th, 2017)
			// The calculation is done by direct integral using MATLAB and fit the result into polynomial (5th order). The variety range is:
			// Wtd = 0.026 to 1 eV.        E_peakp-phip = -1 to 1 eV.     E_peakp = Epeak-Ev
			double Ep = phip-E_peakp;   // Ep = Efp-Epeak
			if (fabs(Ep)>1)
				cout<<"WARNING: The difference between E_peakp and Efp is too large, which may cause inaccuracy."<<endl;
			double p_trap = - 0.0002766 - 0.08858*Ep + 1*Wtd - 0.004624*pow(Ep,2) - 2.973*Ep*Wtd - 0.00552*pow(Wtd,2) + 0.4243*pow(Ep,3) + 0.0116*pow(Ep,2)*Wtd + 2.744*Ep*pow(Wtd,2)
					- 0.03929*pow(Wtd,3) + 0.007633*pow(Ep,4) + 2.475*pow(Ep,3)*Wtd + 0.05881*pow(Ep,2)*pow(Wtd,2) + 4.609*Ep*pow(Wtd,3) + 0.1758*pow(Wtd,4)
					- 0.3877*pow(Ep,5) - 0.03383*pow(Ep,4)*Wtd - 4.877*pow(Ep,3)*pow(Wtd,2) - 0.05436*pow(Ep,2)*pow(Wtd,3) - 5.09*Ep*pow(Wtd,4) - 0.214*pow(Wtd,5);
        		holeDensity = 2*gvV*pow( (mvEff*Vth*ChargeQ/(2*Pi*SQUARE(Planckba))), 1.5) * fermiIntegralHalf(-phip/Vth) + Ntd * p_trap;
				// + Ntd*exp(-(phip-E_peakp)/(nwp*Vth))*(nwp*Vth)*ChargeQ*log1p(exp((phip-E_peakp)/(nwp*Vth)));
		}
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
// Shouldn't include Ntd and Nta in charge density calculation
    return ( Nd - Na  + holeDensity(phip) - electronDensity(phin) );
};

/**
 * return the total (fixed) charge density for insulator
 */
double Material::fixChargeDensity()
{
	if (type != Dielectric) {
		std::cerr << " Error: fixChargeDensity() method is only defined for insulator" << std::endl;
	}
// Shouldn't include Ntd and Nta in fix charge density calculation
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

