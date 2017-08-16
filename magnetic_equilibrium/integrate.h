#ifndef GUARD_integrate_h
#define GUARD_integrate_h


#include <vector>


template <typename Deriv_object>
void integrate(std::vector<std::vector<double> >& r,
	std::vector<std::vector<double> >& dr,
	Deriv_object& deriv, double ti, double tf, double dt)
{
	while( (ti+dt) <= tf){
		deriv(r,dr,dt);
		ti += dt;
	}
	if(ti < tf) {
		deriv(r,dr,tf-ti);
		ti += dt;
	}
}


template <typename Deriv_object>
void integrate_nof(std::vector<std::vector<double> >& r,
	std::vector<std::vector<double> >& dr,
	Deriv_object& deriv, double ti, double tf, double dt)
{
	while( (ti+dt) <= tf){
		deriv.nof(r,dr,dt);
		ti += dt;
	}
	if(ti < tf) {
		deriv.nof(r,dr,tf-ti);
		ti += dt;
	}
}

/*
template <typename Deriv_object>
void integrate(std::vector<std::vector<double> >& r,
	std::vector<std::vector<double> >& dr,
	Deriv_object& deriv, double ti, double tf, double dtmax)
{
	double dt = dtmax;
	while( (ti+dt) <= tf){
		deriv(r,dr,dt);
		ti += dt;
	}
	if(ti < tf) {
		integrate(r,dr,deriv,ti,tf,tf-ti);
	}
}
*/
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
