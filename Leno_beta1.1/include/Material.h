/*
 * Material.h
 *
 *  Created on: Oct 8, 2015
 *      Author: oscar
 */

#ifndef MATERIAL_H_
#define MATERIAL_H_

#include "Params.h"

class Material : public Params {

/** Set to public because material params are unique
 *  in this program.
 */
public:
    int index;
    int type;
    std::string name;
    double dieleY;
    double dieleX;
    double electronAffinity; // postive number
    double bandGap;
    double Nd;
    double Na;
    int dimension;
    double t2D; // the thickness of 2D layer (ONLY for 2D)
    double mcEff;
    double mvEff;
    int gvC; // vally degeneracy
    int gvV;

public:
    Material(); // // Every constructor in the inheritance hierarchy gets called, in the order Base -> Derived. Destructors get called in the reverse order.
    Material(double _T, int _index, int _type, const char *pname, double _dielectricConstant, double _dieleIP, double _electronAffinity, double _bandGap, int _dimension, \
        	double _t2D, double _Nd, double _Na, double _mcEff, double _mvEff, int _gvC, int _gvV);
    Material(int _index, int _type, const char *pname, double _dielectricConstant, double _dieleIP, double _electronAffinity, double _bandGap, int _dimension, \
    		double _t2D, double _Nd, double _Na, double _mcEff, double _mvEff, int _gvC, int _gvV);
    virtual ~Material();

public:
	double electronDensity(double phin);
	double holeDensity(double phip);
    double chargeDensity(double phin, double phip);
    double quantumCapa(double phin, double phip);

};

#endif /* MATERIAL_H_ */
