/*
 * Transport.h
 *
 *  Created on: Oct 12, 2015
 *      Author: oscar
 */

#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include "Poisson1D.h"
class Transport : public Params {
private:
	std::vector<double> interTunnelCurrent, likeTunnelCurrent;
	double 	mb0Square, correlationLengthSquare, broadenConstantSquare, ovSquare;

public:
	Transport();
	virtual ~Transport();

	double interTunnel2DOFP(double Ect, double Evt, double Ecb, double Evb, double Vds, double precision); // keep constant with origin version
	double likeTunnel2DOFP(double Ect, double Evt, double Ecb, double Evb, double Vds, double precision);

	double hpDDCurrent(double potential_i, double potential_iplus1, double carrierDensity_i, double carrierDensity_iplus1);

	std::vector<double> getInterTunnelCurrent();
	std::vector<double> getLikeTunnelCurrent();

private:
	/*
	 * static member function is int (*)(...)
	 * while non-static member funciton is int (Transport::*)(...)
	 */
	// static int fp(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
	// static int f_likeBand(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
};

#endif /* TRANSPORT_H_ */
