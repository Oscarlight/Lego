/*
 * ConstParam.h
 *
 *  Created on: Oct 8, 2015
 *      Author: oscar
 */

#ifndef PARAMS_H_
#define PARAMS_H_

// included dependencies
#include "MathFunct.h";

class Params : public MathFunct {

public:
	// Natural unit:
	const double nat_T = 300; // default Temp.
	const double nat_m0 = 9.1e-31; // kg
	const double nat_Planck = 6.62606957e-34; // J*s
	const double nat_Planckba = Planck/(2*Pi);
	const double nat_Planck_ev = 4.135667517e-15; // eV*s
	const double nat_ChargeQ = 1.6e-19; // C
	const double nat_E0 = 8.854187817e-12; // F/m
	const double nat_kb = 1.3806488e-23; // m2*kg*s-2*K-1

protected:
	// Lego unit:
	const double Pi = 3.1415926;
	const double kb;
	const double T;
	const double Vth;
	const double m0;
	const double Planck;
	const double Planckba;
	const double Planck_ev;
	const double Planckba_ev;
	const double ChargeQ;
	const double E0;

protected:
	// System const:
	const double LARGE = 1E100; // a large number 1E100
	const int    Dirichlet=1;
	const int    Neumann=2;
	const int    Semiconductor=3;
	const int    Dielectric=4;


public:
	Params();
	Params(double T);
	virtual ~Params();

protected:
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
