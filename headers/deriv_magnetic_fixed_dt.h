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
			double vv, double qBB, double omegaa,
			double betaa, double epss, double sigmaa,
			int seedd):
			N(NN),L(LL), Dt(Dtt), sqrt_2Dt(std::sqrt(2*Dtt)),
			sqrt_2Dr(std::sqrt(2*Drr)),
			v(vv),qB(qBB),omega(omegaa),
			beta(betaa), eps(epss), sigma(sigmaa),
			sigma6(pow(sigma,6.)),seed(seedd), generator(seed),ndist(0.,1.),
			F(N,std::vector<double>(3,0.))
			{}


	// calculate derivatives and new positions
	void operator() (
			std::vector<std::vector<double> >& r,
			std::vector<std::vector<double> >& dr,
			std::vector<std::vector<double> >& p,
			std::vector<std::vector<double> >& dp,
			double dt,bool err, double maxForce);

	int get_N() { return N;}
	double get_L() { return L;}
	double get_Dt() { return Dt;}
	double get_Dr() { return 0.5*sqrt_2Dr*sqrt_2Dr;}
	double get_v() { return v;}
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
	double v;
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
	void update_F(const std::vector<std::vector<double> >& r);
	// force between two particles
	double f(const double&);

	//position dependent magnetic field
	double wc(const std::vector<double>& ri)
		{return qB*std::sin(omega*ri[1]);}
	double wcp(const std::vector<double>& ri)
		{return qB*omega*std::cos(omega*ri[1]);}

	// random number generator
	std::default_random_engine generator;
	std::normal_distribution<double> ndist;
};



double Deriv::f(const double& r)
{
	double sr6 = sigma6/(r*r*r*r*r*r);
	return eps*(48.*sr6*sr6-24*sr6)/(r*beta);
}

void Deriv::update_F(
	const std::vector<std::vector<double> >& r)
{
	if( eps > 0 ) {
		double abs_r,abs_f,dx,dy,dz;

		// set force on particle i to zero
		for(int i=0;i<N;++i)
			std::fill(F[i].begin(),F[i].end(),0.);

		for(int i=0;i<N;++i) {

			// add force on i due to j
			for(int j=i+1;j<N;++j) {
				dx = r[j][0] - r[i][0];
				dy = r[j][1] - r[i][1];
				dz = r[j][2] - r[i][2];
				dx -= L*round(dx/L);
				dy -= L*round(dy/L);
				dz -= L*round(dz/L);
				abs_r = sqrt(dx*dx+dy*dy+dz*dz);
				if(abs_r < sigma*pow(2.,1./6) ) {
					abs_f = f(abs_r)/abs_r;
					F[i][0] -= abs_f*dx;
					F[i][1] -= abs_f*dy;
					F[i][2] -= abs_f*dz;
					F[j][0] += abs_f*dx;
					F[j][1] += abs_f*dy;
					F[j][2] += abs_f*dz;
				}
			}

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
		double dt,bool err, double maxForce)
{

	double sqrt_dt = std::sqrt(dt);
	double etaX, etaY, etaZ;	// random numbers for the orientation vector
	double wci;	// position dep. wc
	double wcip;
	double D;
	double drx, dry;

	
	for(int i=0;i<N;++i) {

		//check if forces do not exceed critical value			
		//assert(abs(F[i][0])*dt<sigma);
		//assert(abs(F[i][1])*dt<sigma);
		//assert(abs(F[i][2])*dt<sigma);

		
		if( v > 0 ) {
			// calculate r increment

			if(qB>0) {
				dr[i][0] = (v*p[i][0] + F[i][0])*dt + ndist(generator)*sqrt_dt*sqrt_2Dt;
				dr[i][1] = (v*p[i][1] + F[i][1])*dt + ndist(generator)*sqrt_dt*sqrt_2Dt;
				dr[i][2] = (v*p[i][2] + F[i][2])*dt + ndist(generator)*sqrt_dt*sqrt_2Dt;
				wci = wc(r[i]);
				wcip = wcp(r[i]);
				D = 1. + wci*wci;
				drx = dr[i][0];
				dry = dr[i][1];
				double len = drx*drx+dry*dry;
				// act with Q
				dr[i][0] = (1-wci*wci/D)*drx - (wci/D)*dry;
				dr[i][1] = (1-wci*wci/D)*dry + (wci/D)*drx;
				double len2 = dr[i][0]*dr[i][0]+dr[i][1]*dr[i][1];
				dr[i][0] *= std::sqrt(len/len2);
				dr[i][1] *= std::sqrt(len/len2);
				// add A (preserves eq. dist.)
				dr[i][0] -= wcip*(1-wci*wci)*dt/(D*D);
				dr[i][1] -= 2.*wci*wcip*dt/(D*D);
			} else {
				dr[i][0] = (v*p[i][0] + F[i][0])*dt + ndist(generator)*sqrt_dt*sqrt_2Dt;
				dr[i][1] = (v*p[i][1] + F[i][1])*dt + ndist(generator)*sqrt_dt*sqrt_2Dt;
				dr[i][2] = (v*p[i][2] + F[i][2])*dt + ndist(generator)*sqrt_dt*sqrt_2Dt;
			}
			
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

		} else {

			if(qB>0) {
				dr[i][0] = ndist(generator)*sqrt_dt*sqrt_2Dt + F[i][0]*dt;
				dr[i][1] = ndist(generator)*sqrt_dt*sqrt_2Dt + F[i][1]*dt;
				dr[i][2] = ndist(generator)*sqrt_dt*sqrt_2Dt + F[i][2]*dt;
				wci = wc(r[i]);
				wcip = wcp(r[i]);
				D = 1. + wci*wci;
				drx = dr[i][0];
				dry = dr[i][1];

				double len = drx*drx+dry*dry;
				// act with Q
				dr[i][0] = (1-wci*wci/D)*drx - (wci/D)*dry;
				dr[i][1] = (1-wci*wci/D)*dry + (wci/D)*drx;
				double len2 = dr[i][0]*dr[i][0]+dr[i][1]*dr[i][1];
				dr[i][0] *= len/len2;
				dr[i][1] *= len/len2;
				// add A (preserves eq. dist.)
				dr[i][0] -= wcip*(1-wci*wci)*dt/(D*D);
				dr[i][1] -= 2.*wci*wcip*dt/(D*D);				
			} else {
				dr[i][0] = ndist(generator)*sqrt_dt*sqrt_2Dt + F[i][0]*dt;
				dr[i][1] = ndist(generator)*sqrt_dt*sqrt_2Dt + F[i][1]*dt;
				dr[i][2] = ndist(generator)*sqrt_dt*sqrt_2Dt + F[i][2]*dt;
			}
			r[i][0] += dr[i][0];
			r[i][1] += dr[i][1];
			r[i][2] += dr[i][2];
		}

	}
}


#endif
