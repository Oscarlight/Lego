/*
 * Tunnelling.cpp
 *
 *  Created on: Oct 17, 2015
 *      Author: oscar
 */

#include "Tunnelling.h"

Tunnelling::Tunnelling() {
	mb0Square = SQUARE(0.02);
	correlationLengthSquare = SQUARE(10*cmNL(1E-7));
	broadenConstantSquare = SQUARE(0.01);
    ovSquare = 0.0278;
}

Tunnelling::~Tunnelling() {
	// TODO Auto-generated destructor stub
}

// Tunnel current
// p-type
/*
 * Don't add static keyword here -> cannot declare member function to have static linkage
 * static function cannot access non-static members.
 */
int f_interband(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
	//
	double * currentParaArray = (double *) fdata;
	// --------- name list ------------
    double mb0Square = currentParaArray[0]; // Vrms^2
    double ovSquare = currentParaArray[1]; // exp(-decay*di)^2
    double correlationLengthSquare = currentParaArray[2];
    double broadenConstantSquare = currentParaArray[3];
    double Ect = currentParaArray[4];
    double Evt = currentParaArray[5];
    double Ecb = currentParaArray[6];
    double Evb = currentParaArray[7];
    double mvt = currentParaArray[8];
    double mcb = currentParaArray[9];
    double mct = currentParaArray[10];
    double mvb = currentParaArray[11];
    double Vds = currentParaArray[12];
    double mink = currentParaArray[13];

	double Planckba = currentParaArray[14];
	double Planckba_ev = currentParaArray[15];
	double ChargeQ = currentParaArray[16];
	double Pi = currentParaArray[17];
	double Vth = currentParaArray[18];


	double kt = x[0];
	double kb = x[1];
	double dtheta = x[2];
	// -------- end name list ---------

	double Eckf, Evkf, dEkf;
	double Ef_f, Ef_i;
	if ( std::abs(Ect - Evb)  >  std::abs(Evt - Ecb)) {  // p type
		Eckf = Ecb + SQUARE(Planckba) * SQUARE(kt)/(2*mct*ChargeQ);
		Evkf = Evt - SQUARE(Planckba) * SQUARE(kb)/(2*mvb*ChargeQ);
		Ef_f = Eckf;
		Ef_i = Evkf + Vds;
	} else { // n type
		Eckf = Ect + SQUARE(Planckba) * SQUARE(kt)/(2*mct*ChargeQ);
		Evkf = Evb - SQUARE(Planckba) * SQUARE(kb)/(2*mvb*ChargeQ);
		Ef_f = Eckf + Vds;
		Ef_i = Evkf;
	}
	dEkf = Eckf - Evkf;
	// DOS Broadening
	double Pq = (1/std::sqrt(2*Pi*broadenConstantSquare)) * std::exp(-SQUARE(dEkf)/(2*broadenConstantSquare));
	// Fermi Distribution
	double FDD = 1/(1 + std::exp(Ef_i/Vth)) - 1/(1 + std::exp(Ef_f/Vth));
	double FDD0 = 1/(1 + std::exp(Evkf/Vth)) - 1/(1 + std::exp(Eckf/Vth));
	FDD -= FDD0;

	// Scattering
	double Sq = Pi * mb0Square * correlationLengthSquare/std::pow( (1 + 0.5 * correlationLengthSquare * (SQUARE(kt) + SQUARE(kb) + 2*kt*kb*std::cos(dtheta))), 1.5);
	// current
	double gv = 2;
	fval[0] = gv * ((4 * Pi * ChargeQ)/Planckba_ev) * ovSquare * ( 1/std::pow((2*Pi), 4) ) * 2 * Pi * Sq * Pq * FDD * kt * kb; // A/m^2
	return (0); // success
}

int f_likeBand(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) // didn't modify yet
{
    //
    double * currentParaArray = (double *) fdata; // maybe wrong
    // --------- name list ------------
    double mb0Square = currentParaArray[0]; // Vrms^2
    double ovSquare = currentParaArray[1]; // exp(-decay*di)^2
    double correlationLengthSquare = currentParaArray[2]; // 1/m^2
    double broadenConstantSquare = currentParaArray[3]; // eV
    double Ect = currentParaArray[4];
    double Evt = currentParaArray[5];
    double Ecb = currentParaArray[6];
    double Evb = currentParaArray[7];
    double mvt = currentParaArray[8];
    double mcb = currentParaArray[9];
    double mct = currentParaArray[10];
    double mvb = currentParaArray[11];
    double Vds = currentParaArray[12];
    double mink = currentParaArray[13];

	double Planckba = currentParaArray[14];
	double Planckba_ev = currentParaArray[15];
	double ChargeQ = currentParaArray[16];
	double Pi = currentParaArray[17];
	double Vth = currentParaArray[18];

    double kt = x[0];
    double kb = x[1];
    double dtheta = x[2];
    // -------- end name list ---------
    // conduction-conduction band
    double Eckb = Ecb + SQUARE(Planckba) * SQUARE(kb)/(2*mcb*ChargeQ);
    double Eckt = Ect + SQUARE(Planckba) * SQUARE(kt)/(2*mct*ChargeQ);
    double dEck = Eckb - Eckt;
    // valence-valence band
    double Evkb = Evb - SQUARE(Planckba) * SQUARE(kb)/(2*mvb*ChargeQ);
    double Evkt = Evt - SQUARE(Planckba) * SQUARE(kt)/(2*mvt*ChargeQ);
    double dEvk = Eckb - Eckt;

    // DOS Broadening
    double Pq_c = (1/std::sqrt(2*Pi*broadenConstantSquare)) * std::exp(-SQUARE(dEck)/(2*broadenConstantSquare));
    double Pq_v = (1/std::sqrt(2*Pi*broadenConstantSquare)) * std::exp(-SQUARE(dEvk)/(2*broadenConstantSquare));
    // testing
    // double Pq_c = 1, Pq_v = 1;

    // Fermi Distribution
    double Ecfb = Eckb;
    double Ecft = Eckt + Vds;
    double FDD_c = 1/(1 + std::exp(Ecft/Vth)) - 1/(1 + std::exp(Ecfb/Vth)); // Boltzman Distribution may work as well
    double FDD_c0 = 1/(1 + std::exp(Eckt/Vth)) - 1/(1 + std::exp(Eckb/Vth));
    FDD_c -= FDD_c0;

    double Evfb = Evkb;
    double Evft = Evkt + Vds;
    double FDD_v = 1/(1 + std::exp(Evft/Vth)) - 1/(1 + std::exp(Evfb/Vth)); // how to reverse +/- for hole branch current, may be wrong
    double FDD_v0 = 1/(1 + std::exp(Evkt/Vth)) - 1/(1 + std::exp(Evkb/Vth));
    FDD_v -= FDD_v0;

    // Scattering
    // double Sq_c = Pi * mb0Square * correlationLengthSquare/std::pow( (1 + 0.5 * correlationLengthSquare * (SQUARE(kt-mink) + SQUARE(kb) + 2*(kt-mink)*kb*std::cos(dtheta))), 1.5);
    // double Sq_v = Pi * mb0Square * correlationLengthSquare/std::pow( (1 + 0.5 * correlationLengthSquare * (SQUARE(kt) + SQUARE(kb-mink) + 2*kt*(kb-mink)*std::cos(dtheta))), 1.5);
    double Sq = Pi * mb0Square * correlationLengthSquare/std::pow( (1 + 0.5 * correlationLengthSquare * (SQUARE(kt) + SQUARE(kb) + 2 * kt * kb * std::cos(dtheta))), 1.5);

    // current
    double gv = 2;
    fval[0] = gv * ((4 * Pi * ChargeQ)/Planckba_ev) * ovSquare * ( 1/std::pow((2*Pi), 4) ) * 2 * Pi * Sq * (Pq_c * FDD_c + Pq_v * FDD_v) * kt * kb; // A/m^2
    // CAUTION: the sign in FDD_v and the sign within Pq_c * FDD_c + Pq_v * FDD_v may be WRONG, keep eye on it!

    return (0); // success
}


//// half-point DD current (Scharfetter-Gummel Method)
//double Transport::hpDDCurrent(double potential_i, double potential_iplus1, double carrierDensity_i, double carrierDensity_iplus1)
//{
//	double hpJ;
//	double mobility = 0.1; // m2/Vs * 10000 = cm2/Vs
//	double t = (potential_iplus1 - potential_i)/Vth;
//	/* CAUTION! Eintein Relation is used w/o verificaiton */
//	hpJ = ChargeQ * mobility * Vth / (meshx * A_TO_METER) * ( bernoulli(t) * carrierDensity_iplus1 - bernoulli(-t) * carrierDensity_i ); // A/m
//	return hpJ;
//}

// May be a little unnecessary complicated arguments, in consistent with orignal matlab code
// largely copied from the original Thin-TFET C++ code
double Tunnelling::interTunnel2DOFP(double Ect, double Evt, double Ecb, double Evb, double Vds, double precision)
{

    double mvt = 0.4*m0, mcb = 0.3*m0;
    double mct = 0.3*m0, mvb = 0.4*m0;

    // self-adapting range for maxk, (if needed, can implement the same method on mink)
    double dEcv = std::min(std::abs(Evb - Ect), std::abs(Evt - Ecb)) + 0.2; // 0.2 is the margin, for p-type Thin-TFET, 0.2 -> 0.3 no different, 0.2 -> 0.1 < 1% difference
    double maxk = std::sqrt( 2 * mct * dEcv * ChargeQ)/Planckba;

    // parameter
    double mink = 0; // look at the thermal tail part
	double xmin[3] = {mink,mink,0}, xmax[3] = {maxk,maxk,2*Pi}, val, err;

	// for p-type TFET
//	double qphip_t = Eg_t - qphin_t;
	double fdata[] = {mb0Square, ovSquare, correlationLengthSquare, broadenConstantSquare, Ect, Evt, Ecb, Evb, mvt, mcb, mct, mvb, Vds, mink, Planckba, Planckba_ev, ChargeQ, Pi, Vth};
    pcubature(1, f_interband, fdata, 3, xmin, xmax, 0, 0, precision, ERROR_INDIVIDUAL, &val, &err);
    //

    printf("inter-band tunnel current = %0.10g +/- %g A/m \n", val * cLN(1)/sLN(1) / SQUARE(cmLN(1)) * 1E4 * 15E-9,
    		err * cLN(1)/sLN(1) / SQUARE(cmLN(1)) * 1E4 * 15E-9);

    interTunnelCurrent.push_back(val);

    return (val);
}




double Tunnelling::likeTunnel2DOFP(double Ect, double Evt, double Ecb, double Evb, double Vds, double precision)
{
    double mvt = 0.4*m0, mcb = 0.3*m0;
    double mct = 0.3*m0, mvb = 0.4*m0;

    // self-adapting range for maxk, (if needed, can implement the same method on mink)
    double dEc = std::abs(Ect - Ecb) + 0.2; // 0.2 is the margin
    double maxk_c = std::sqrt( 2 * mct * dEc * ChargeQ)/Planckba;
    double dEv = std::abs(Evt - Evb) + 0.2;
    double maxk_v = std::sqrt( 2 * mvt * dEv * ChargeQ)/Planckba;

    // parameter
    double mink = 0; // look at the thermal tail part
    double maxk = std::max(maxk_v, maxk_c);
    // CAUSION: if maxk is too large, the integral will become 0 OR sometime the error will increase (i.e. abrupt change)
    // If it is too small, it won't cover the interesting range
    // testing is necessary to find the opimal maxk
    double xmin[3] = {mink,mink,0}, xmax[3] = {maxk,maxk,2*Pi}, val, err;

    double fdata[] = {mb0Square, ovSquare, correlationLengthSquare, broadenConstantSquare, Ect, Evt, Ecb, Evb, mvt, mcb, mct, mvb, Vds, mink, Planckba, Planckba_ev, ChargeQ, Pi, Vth};
    pcubature(1, f_likeBand, fdata, 3, xmin, xmax, 0, 0, precision, ERROR_INDIVIDUAL, &val, &err);
    // accuracy also important
    //
    printf("like-band tunnel current = %0.10g +/- %g A/m \n", val * cLN(1)/sLN(1) / SQUARE(cmLN(1)) * 1E4 * 15E-9,
    		err * cLN(1)/sLN(1) / SQUARE(cmLN(1)) * 1E4 * 15E-9 );

    likeTunnelCurrent.push_back(val);

    return (val);
}


std::vector<double> Tunnelling::getInterTunnelCurrent() {
		return (interTunnelCurrent);
}


std::vector<double> Tunnelling::getLikeTunnelCurrent() {
		return (likeTunnelCurrent);
}

