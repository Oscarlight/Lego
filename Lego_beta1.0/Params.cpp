/*
 * ConstParam.cpp
 *
 *  Created on: Oct 8, 2015
 *      Author: oscar
 */

#include "Params.h"

Params::Params() {
	std::cout << "Default constructor of Params is called with T = 300 K" << std::endl;
	init();
}

Params::Params(double _T) {
	// T: temperature
	std::cout << "constructor of Params is called to set T" << std::endl;
	nat_T = _T;
	init(); 
}

Params::~Params() {}

/**
 * set (1E18 cm^-3)^(-1/3) as the Lego unit of cm
 */
double Params::cmNL(double a) {
	return (a * 1E6);
}

double Params::cmLN(double a) {
	return (a * 1E-6);
}

/**
 * set 10^12 s as Lego unit of s
 */
double Params::sNL(double a) {
	return (a * 1E12);
}

double Params::sLN(double a) {
	return (a * 1E-12);
}

/**
 * set 1/1.6E-19 C as unit of charge
 */
double Params::cNL(double a) {
	return (a / nat_ChargeQ);
}

double Params::cLN(double a) {
	return (a * nat_ChargeQ);
}

/**
 * set 1/300 K as the unit for temperature
 */
double Params::tNL(double a) {
	return (a / 300);
}

double Params::tLN(double a) {
	return (a * 300);
}

/**
 * J = C*V
 */
double Params::jNL(double a) {
	return (a * cNL(1) * 1);
}

/**
 * kg = J/(m^2*s^2)
 */
double Params::kgNL(double a) {
	return (a * jNL(1)/ SQUARE(cmNL(100) * sNL(1)));
}

/**
 * init: the Natural unit to convert to Lego unit
 */
void Params::init() {
	// m2*kg*s-2*K-1
	kb = nat_kb * cmNL(100) * kgNL(1) / ( SQUARE(sNL(1)) * tNL(1) );
	T = tNL(nat_T);
	// Vth = kb T / q
	Vth = kb * T / (cNL(1));
	// m0
	m0 = nat_m0 * kgNL(1);
	// m^2*kg/s == J*s
	Planck = nat_Planck * jNL(1) * sNL(1);
	Planckba = Planck/(2*Pi);
	// eV*s
	Planck_ev = nat_Planck_ev * 1 * sNL(1);
	Planckba_ev = Planck_ev/(2*Pi);
	// electron charge: e
	ChargeQ = nat_ChargeQ * cNL(1);
	// F/m, F = C/V
	E0 = nat_E0 * cNL(1) / ( 1 * cmNL(100) );
}
