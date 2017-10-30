#ifndef GUARD_integrate_h
#define GUARD_integrate_h


#include <vector>

template <typename Deriv_object>
void integrate(std::vector<std::vector<double> >& r,
	std::vector<std::vector<double> >& dr,
	std::vector<std::vector<double> >& p,
	std::vector<std::vector<double> >& dp,
	Deriv_object& deriv, double ti, double tf, double dt)
{
	while( (ti+dt) <= tf){
		deriv(r,dr,p,dp,dt);
		ti += dt;
	}
	if(ti < tf) {
		integrate(r,dr,p,dp,deriv,ti,tf,tf-ti);
	}
}

template <typename Deriv_object>
void integrate(std::vector<std::vector<double> >& r,
	std::vector<std::vector<double> >& dr,
	std::vector<std::vector<double> >& v,
	std::vector<std::vector<double> >& dv,
	std::vector<std::vector<double> >& p,
	std::vector<std::vector<double> >& dp,
	Deriv_object& deriv, double ti, double tf, double dt)
{
	while( (ti+dt) <= tf){
		deriv(r,dr,v,dv,p,dp,dt);
		ti += dt;
	}
	if(ti < tf) {
		integrate(r,dr,v,dv,p,dp,deriv,ti,tf,tf-ti);
	}
}

#endif
