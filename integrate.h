#ifndef GUARD_integrate_h
#define GUARD_integrate_h


#include <vector>


#include "deriv.h"


void integrate(std::vector<std::vector<double> >& r,
	std::vector<std::vector<double> >& dr,
	std::vector<std::vector<double> >& p,
	std::vector<std::vector<double> >& dp,
	Deriv& deriv, double ti, double tf, double dt)
{
	while( (ti+dt) < tf){
		deriv(r,dr,p,dp,dt);
		ti += dt;
	}
	if(ti < tf)
		deriv(r,dr,p,dp,tf-ti);
}


void integrate(std::vector<std::vector<double> >& r,
	std::vector<std::vector<double> >& p,
	Deriv& deriv, double ti, double tf, double dt)
{
	std::vector<std::vector<double> > dr = r;
	std::vector<std::vector<double> > dp = p;
	while( (ti+dt) < tf){
		deriv(r,dr,p,dp,dt);
		ti += dt;
	}
	if(ti < tf)
		deriv(r,dr,p,dp,tf-ti);
}




#endif
