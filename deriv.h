#ifndef GUARD_deriv_h
#define GUARD_deriv_h

/* 
 * object with the SDE
 * the integrarion within a sphere
 * is done in integration.h
 */


#include <cmath>
#include <vector>
#include <random>

#include "vecmanip.h"

struct Deriv {
	public:

	Deriv(int NN, double Dtt, double Drr,double betaa,
			double epss, double sigmaa, int seedd):
			N(NN),sqrt_2Dt(std::sqrt(2*Dtt)), sqrt_2Dr(std::sqrt(2*Drr)),
			beta(betaa), eps(epss), sigma(sigmaa),
			seed(seedd), generator(seed),ndist(0.,1.){}


	// calculate derivatives and new positions
	void operator() (
			std::vector<std::vector<double> >& r,
			std::vector<std::vector<double> >& dr,
			std::vector<std::vector<double> >& p,
			std::vector<std::vector<double> >& dp,
			double dt);

	int get_N() { return N;}
	double get_Dt() { return 0.5*sqrt_2Dt*sqrt_2Dt;}
	double get_Dr() { return 0.5*sqrt_2Dr*sqrt_2Dr;}
	double get_beta() { return beta;}
	double get_eps() { return eps;}
	double get_sigma() {return sigma;}
	std::vector<std::vector<std::vector<double> > > get_F() {return F;}

	private:
	int N;				// number of particles
	double sqrt_2Dt;	// sqrt(2*Dt)
	double sqrt_2Dr;	// sqrt(2*Dr)
	double beta;
	double eps;
	double sigma;
	int seed;			// seed for the random generator

	// the force matrix
	std::vector<std::vector<std::vector<double> > > F;
	// calculate the force matrix
	void update_F(const std::vector<std::vector<double> >& r);
	// force between two particles
	std::vector<double> f(const std::vector<double>& r1,
			const std::vector<double>& r2);

	// random number generator
	std::default_random_engine generator;
	std::normal_distribution<double> ndist;
};

std::vector<double> Deriv::f(const std::vector<double>& r1,
		const std::vector<double>& r2)
{
	std::vector<double> r = subtr(r1,r2);
	double abs_r = len_vec(r);
	double abs_f = 0;
	if(abs_r<sigma*pow(2,1./6))
		abs_f = 4*eps*(pow(sigma/abs_r,12)-pow(sigma/abs_r,6))/beta;
	mult_by(r,abs_f/sqrt(abs_r));

	return r;
}

void Deriv::update_F(
		const std::vector<std::vector<double> >& r)
{
	for(int i=0;i<N;++i) {
		for(int j=i+1;j<N;++j) {
			F[i][j] = f(r[i],r[j]); 
			F[j][i] = mult(F[i][j],-1.);
		}
	}
}

// The () operator calculates the increment in r and p (dr and dp) at r,p
// for a time increment dt and adds it to r and p
void Deriv::operator() (
		std::vector<std::vector<double> >& r,
		std::vector<std::vector<double> >& dr,
		std::vector<std::vector<double> >& p,
		std::vector<std::vector<double> >& dp,
		double dt)
{
	double sqrt_dt = std::sqrt(dt);
	double etaX, etaY, etaZ;
	update_F(r);
	for(int i=0;i<N;++i) {
		// calculate r increment
		dr[i][0] = ndist(generator)*sqrt_dt*sqrt_2Dt;
		dr[i][1] = ndist(generator)*sqrt_dt*sqrt_2Dt;
		dr[i][2] = ndist(generator)*sqrt_dt*sqrt_2Dt;
		// add inter particle forces
		for(int j=0;j<N;++j) { // loop over j=i;j<N;++j ?
			if(i==j) continue;
			dr[i][0] += F[i][j][0];
			dr[i][1] += F[i][j][1];
			dr[i][2] += F[i][j][2];
		}

		// calculate p increment
		etaX = ndist(generator)*sqrt_dt*sqrt_2Dr;
		etaY = ndist(generator)*sqrt_dt*sqrt_2Dr;
		etaZ = ndist(generator)*sqrt_dt*sqrt_2Dr;
		dp[i][0] = (etaY*p[i][2] - etaZ*p[i][1]);
		dp[i][1] = (etaZ*p[i][0] - etaX*p[i][2]);
		dp[i][2] = (etaX*p[i][1] - etaY*p[i][0]);
	}
		
}

#endif
