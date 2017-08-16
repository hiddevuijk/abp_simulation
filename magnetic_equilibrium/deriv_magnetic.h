#ifndef GUARD_deriv_magnetic_h
#define GUARD_deriv_magnetic_h



#include <cmath>
#include <vector>
#include <random>
#include <assert.h>
#include "vecmanip.h"


struct Deriv {
	public:

	Deriv(int NN, double LL, double Dtt,
			double qBB, double omegaa,
			double gammaa,int seedd):
			N(NN),L(LL), sqrt_2Dt(std::sqrt(2*Dtt)),
			qB(qBB),omega(omegaa),gamma(gammaa),
			seed(seedd), generator(seed),ndist(0.,1.)
			{}


	// calculate derivatives and new positions
	void operator() (
			std::vector<std::vector<double> >& r,
			std::vector<std::vector<double> >& dr,
			double dt);
	void nof (
			std::vector<std::vector<double> >& r,
			std::vector<std::vector<double> >& dr,
			double dt);


	int get_N() { return N;}
	double get_L() { return L;}
	double get_Dt() { return 0.5*sqrt_2Dt*sqrt_2Dt;}
	double get_qB() { return qB;}
	double get_omega(){return omega;}
	double get_gamma() { return gamma;}

	private:
	int N;			// number of particles
	double L;		// size of the box
	double sqrt_2Dt;	// sqrt(2*Dt)
	double qB;
	double omega;
	double gamma;
	int seed;		// seed for the random generator


	//position dependent magnetic field
	double wc(const std::vector<double>& ri, double qB);
	double wcp(const std::vector<double>& ri, double qB);
	void Q(std::vector<double>& dri, double wci);
	void f1(std::vector<double>& dri, double wci, double wcip,double dt);
	// random number generator
	std::default_random_engine generator;
	std::normal_distribution<double> ndist;
};



double Deriv::wc(const std::vector<double>& ri, double qB)
{
	//if(omega <=0 ) return qB;
	return qB*std::sin(omega*ri[1]);
}

double Deriv::wcp(const std::vector<double>& ri, double qB)
{
	//if(omega <=0 ) return 0.;
	return qB*omega*std::cos(omega*ri[1]);
}


// The () operator calculates the increment in r and p (dr and dp) at r,p
// for a time increment dt and adds it to r and p
void Deriv::operator() (
		std::vector<std::vector<double> >& r,
		std::vector<std::vector<double> >& dr,
		double dt)
{
	double sqrt_dt = std::sqrt(dt);
	double wci;	// position dep. wc
	double wcip;
	for(int i=0;i<N;++i) {
		dr[i][0] = gamma*ndist(generator)*sqrt_dt*sqrt_2Dt;
		dr[i][1] = gamma*ndist(generator)*sqrt_dt*sqrt_2Dt;
		dr[i][2] = gamma*ndist(generator)*sqrt_dt*sqrt_2Dt;
		wci = wc(r[i],qB);
		Q(dr[i],wci);
		wcip = wcp(r[i],qB);
		f1(dr[i],wci,wcip,dt);
		r[i][0] += dr[i][0];
		r[i][1] += dr[i][1];
		r[i][2] += dr[i][2];

	}
}
void Deriv::nof (
		std::vector<std::vector<double> >& r,
		std::vector<std::vector<double> >& dr,
		double dt)
{
	double sqrt_dt = std::sqrt(dt);
	for(int i=0;i<N;++i) {
		dr[i][0] = gamma*ndist(generator)*sqrt_dt*sqrt_2Dt;
		dr[i][1] = gamma*ndist(generator)*sqrt_dt*sqrt_2Dt;
		dr[i][2] = gamma*ndist(generator)*sqrt_dt*sqrt_2Dt;
		r[i][0] += dr[i][0];
		r[i][1] += dr[i][1];
		r[i][2] += dr[i][2];

	}
}


void Deriv::Q(std::vector<double>& dri, double wci)
{

	double drx = dri[0];
	double dry = dri[1];

	double A = wci/(gamma*gamma+wci*wci);
	double B = 1./gamma-gamma/(gamma*gamma+wci*wci);

	dri[0] = drx/gamma + A*dry - B*drx;
	dri[1] = dry/gamma - A*drx - B*dry; 
	dri[2] /= gamma;

}


void Deriv::f1(std::vector<double>& dr, double wci,double wcip,double dt)
{
	double D = gamma*gamma+wci*wci;
	double a = wci/D;
	double b = gamma/D;

	double A = wcip*(gamma*gamma*-wci*wci)/(D*D);
	double B = -2.*gamma*wcip*wci/(D*D);

	dr[1] += (a*A+b*B)*dt;

}






#endif
