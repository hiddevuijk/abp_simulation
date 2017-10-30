#ifndef GUARD_deriv_magnetic_h
#define GUARD_deriv_magnetic_h



#include <cmath>
#include <vector>
#include <random>
#include <assert.h>
#include "vecmanip.h"


struct Deriv {
	public:

	Deriv(int NN, double LL, double Dtt, double Drr,
			double v00,double mm, double qBB, double omegaa,
			double betaa, double epss, double sigmaa,
			int seedd):
			N(NN),L(LL), Dt(Dtt), sqrt_2Dt(std::sqrt(2*Dtt)),
			sqrt_2Dr(std::sqrt(2*Drr)),
			v0(v00),m(mm),qB(qBB),omega(omegaa),
			beta(betaa), eps(epss), sigma(sigmaa),
			sigma6(pow(sigma,6.)),seed(seedd), generator(seed),ndist(0.,1.),
			F(N,std::vector<double>(3,0.))
			{}


	// calculate derivatives and new positions
	void operator() (
			std::vector<std::vector<double> >& r,
			std::vector<std::vector<double> >& dr,
			std::vector<std::vector<double> >& v,
			std::vector<std::vector<double> >& dv,
			std::vector<std::vector<double> >& p,
			std::vector<std::vector<double> >& dp,
			double dt);

	int get_N() { return N;}
	double get_L() { return L;}
	double get_Dt() { return Dt;}
	double get_Dr() { return 0.5*sqrt_2Dr*sqrt_2Dr;}
	double get_v0() { return v0;}
	double get_qB() { return qB;}
	double get_omega(){return omega;}
	double get_beta() { return beta;}
	double get_eps() { return eps;}
	double get_sigma() {return sigma;}
	std::vector<std::vector<double> > get_F() {return F;}

	private:
	int N;			// number of particles
	double L;		// size of the box
	double Dt;
	double sqrt_2Dt;	// sqrt(2*Dt)
	double sqrt_2Dr;	// sqrt(2*Dr)
	double v0;
	double m;
	double qB;
	double omega;
	double beta;
	double eps;
	double sigma;
	double sigma6;
	int seed;		// seed for the random generator

	// the force matrix
	std::vector<std::vector<double> > F;
	// calculate the force matrix

	//position dependent magnetic field
	double wc(const std::vector<double>& ri)
		{return qB*std::sin(omega*ri[1]);}
	double wcp(const std::vector<double>& ri)
		{return qB*omega*std::cos(omega*ri[1]);}

	// random number generator
	std::default_random_engine generator;
	std::normal_distribution<double> ndist;
};





// The () operator calculates the increment in r and p (dr and dp) at r,p
// for a time increment dt and adds it to r and p
void Deriv::operator() (
		std::vector<std::vector<double> >& r,
		std::vector<std::vector<double> >& dr,
		std::vector<std::vector<double> >& v,
		std::vector<std::vector<double> >& dv,
		std::vector<std::vector<double> >& p,
		std::vector<std::vector<double> >& dp,
		double dt)
{

	double sqrt_dt = std::sqrt(dt);
	double etaX, etaY, etaZ;	// random numbers for the orientation vector
	double wci;	// position dep. wc
	double D;
	double drx, dry;

	
	for(int i=0;i<N;++i) {

		//check if forces do not exceed critical value			
		//assert(abs(F[i][0])*dt<sigma);
		//assert(abs(F[i][1])*dt<sigma);
		//assert(abs(F[i][2])*dt<sigma);
		wci = wc(r[i]);
	
		dr[i][0] = v[i][0]*dt;
		dr[i][1] = v[i][1]*dt;
		dr[i][2] = v[i][2]*dt;

		double vnorm = len_vec(v[i]);
		double vx = v[i][0];
		double vy = v[i][1];
		v[i][0] -= wci*vy*dt/m;
		v[i][1] += wci*vx*dt/m;
		normalize(v[i],vnorm);
		

		dv[i][0] = (-v[i][0]+v0*p[i][0])*dt/m + ndist(generator)*sqrt_dt*sqrt_2Dt/m;
		dv[i][1] = (-v[i][1]+v0*p[i][1])*dt/m + ndist(generator)*sqrt_dt*sqrt_2Dt/m;
		dv[i][2] = (-v[i][2]+v0*p[i][2])*dt/m + ndist(generator)*sqrt_dt*sqrt_2Dt/m;
	

		v[i][0] += dv[i][0];
		v[i][1] += dv[i][1];
		v[i][2] += dv[i][2];
		r[i][0] += dr[i][0];
		r[i][1] += dr[i][1];
		r[i][2] += dr[i][2];

		// calculate p increment
		etaX = ndist(generator)*sqrt_dt*sqrt_2Dr;
		etaY = ndist(generator)*sqrt_dt*sqrt_2Dr;
		etaZ = ndist(generator)*sqrt_dt*sqrt_2Dr;
		p[i][0] += (etaY*p[i][2] - etaZ*p[i][1]);
		p[i][1] += (etaZ*p[i][0] - etaX*p[i][2]);
		p[i][2] += (etaX*p[i][1] - etaY*p[i][0]);

		
		normalize(p[i]);
	}
}


#endif
