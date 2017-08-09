#ifndef GUARD_integrate_h
#define GUARD_integrate_h


#include <vector>

//#include "deriv.h"
template <typename Deriv_object>
void integrate(std::vector<std::vector<double> >& r,
	std::vector<std::vector<double> >& dr,
	std::vector<std::vector<double> >& p,
	std::vector<std::vector<double> >& dp,
	Deriv_object& deriv, double ti, double tf, double dtmax)
{
	bool err = false;	// true if F*dt>sigma
	double maxForce;
	double dt = dtmax;
	while( (ti+dt) <= tf){
		deriv(r,dr,p,dp,dt,err,maxForce);
		if(err) {
			deriv(r,dr,p,dp,dt,err,maxForce);
			dt = 0.5*deriv.get_sigma()/maxForce;
			ti += dt;
		} else {
			ti += dt;
			// double dt, but do not exceed dtmax
			if(dt*2.>dtmax) dt = dtmax;
			else dt *= 2.;
		}
	}
	if(ti < tf) {
		integrate(r,dr,p,dp,deriv,ti,tf,tf-ti);
	}
}

/*
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


*/

#endif
