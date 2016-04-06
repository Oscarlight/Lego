/*
 * ConstParam.h
 *
 *  Created on: Oct 8, 2015
 *      Author: oscar
 */

#ifndef PARAMS_H_
#define PARAMS_H_

// included dependencies
#include "MathFunct.h"

class Params : public MathFunct {

public:
	// Natural unit:
	double Pi = 3.1415926;
	double nat_T = 300; // default Temp.
	double nat_m0 = 9.1e-31; // kg
	double nat_Planck = 6.62606957e-34; // J*s
	double nat_Planckba = nat_Planck/(2*Pi);
	double nat_Planck_ev = 4.135667517e-15; // eV*s
	double nat_Planckba_ev = nat_Planck_ev/(2*Pi);; // eV*s
	double nat_ChargeQ = 1.6e-19; // C
	double nat_E0 = 8.854187817e-12; // F/m
	double nat_kb = 1.3806488e-23; //  J/K

public:
//protected:
	// Lego unit:
	double kb;
	double T;
	double Vth;
	double m0;
	double Planck;
	double Planckba;
	double Planck_ev;
	double Planckba_ev;
	double ChargeQ;
	double E0;

public:
	// System const:
	double LARGE = 1E10; // a large number, used to be too large (1E100) but cause matrix to be singular
	int    Dirichlet=1;
	int    Neumann=2;
	int    Semiconductor=3;
	int    Dielectric=4;


public:
	Params();
	Params(double T);
	virtual ~Params();

//protected:
	// convert Nature unit to Lego unit
	double cmNL(double a); // cm
	double sNL(double a); // second
	double cNL(double a); // charge
	double tNL(double a); // temperature

	// convert Lego unit to Natural unit
	double cmLN(double a); // cm
	double sLN(double a); // second
	double cLN(double a); // charge
	double tLN(double a); // temperature

	// derived unit
	double jNL(double a); // J
	double kgNL(double a); // kg

private:
	// convert all Natural unit to Lego unit
	void init();
};

#endif /* PARAMS_H_ */
