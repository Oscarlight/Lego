/*
 * Tunnelling.h
 *
 *  Created on: Oct 17, 2015
 *      Author: oscar
 */

#ifndef SRC_TUNNELLING_H_
#define SRC_TUNNELLING_H_

#include "Transport.h"
class Tunnelling : public Transport{
private:
	std::vector<double> interTunnelCurrent, likeTunnelCurrent;
	double 	mb0Square, correlationLengthSquare, broadenConstantSquare, ovSquare;

public:
	Tunnelling();
	virtual ~Tunnelling();

	double interTunnel2DOFP(double Ect, double Evt, double Ecb, double Evb, double Vds, double precision); // keep constant with origin version
	double likeTunnel2DOFP(double Ect, double Evt, double Ecb, double Evb, double Vds, double precision);

	double hpDDCurrent(double potential_i, double potential_iplus1, double carrierDensity_i, double carrierDensity_iplus1);

	std::vector<double> getInterTunnelCurrent();
	std::vector<double> getLikeTunnelCurrent();
};

#endif /* SRC_TUNNELLING_H_ */
