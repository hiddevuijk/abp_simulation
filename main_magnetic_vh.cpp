
#include <iostream>
#include <vector>
#include <string>
#include <random>

#include "read.h"
#include "deriv_magnetic.h"
#include "integrate.h"
#include "fcc_lattice.h"
#include "vanHove.h"


using namespace std;


typedef vector<vector<double> > M2;
typedef vector<vector<vector<double> > > M3;

int main(int argc, char *argv[])
{

	int navg;
	int N, seed;
	double Dt,Dr,v0,qB,omega0,beta,eps,sigma,L,dt,tf,teq;
	double rho;
	double bs;		// binsize
	double tbs;		// time-step size
	double tmax;		// end time

	// name of the output file
	string name;

	// name of the input file
	// default is input.txt, otherwise commandline input
	string input_name = "input.txt";
	if(argc == 2) {
		input_name = argv[1];
		input_name += ".txt";
	}

	// read variables from input file
	read_variables_magnetic_vanHove(N,rho,Dt,Dr,v0,qB,omega0,beta,eps,
			sigma,dt,tf,teq,navg,bs,tbs,tmax,seed,name,input_name);

	// L: box size
	L = pow(N/rho,1./3);

	// omega = omega0 2 pi / L
	double omega = omega0*2*acos(-1)/L;

	int Nbin = ceil(L/bs);
	vector<double> bins(Nbin,0.0);
	for(int i=0;i<Nbin;++i)
		bins[i] = (i+.5)*bs;

	int Ntbin = ceil(tmax/tbs);
	vector<double> tbins(Ntbin,0.0);
	for(int i=0;i<Ntbin;++i)
		tbins[i] = (i+.5)*bs;

	M2 r0(N,vector<double>(3));
	M2 r(N,vector<double>(3));
	M2 dr(N,vector<double>(3));
	M2 p(N,vector<double>(3,1.));
	M2 dp(N,vector<double>(3,1.));

	M3 vHx(Ntbin,M2(Nbin,vector<double>(Nbin,0.0)));
	M3 vHy(Ntbin,M2(Nbin,vector<double>(Nbin,0.0)));
	M3 vHz(Ntbin,M2(Nbin,vector<double>(Nbin,0.0)));


	// initalize r vector: put particles on a fcc lattice
	init_r_fcc(r,N,sigma,L);

	// init deriv objec to perform integration
	Deriv deriv(N,L,Dt,Dr,v0,qB,omega,beta,eps,sigma,seed);

	// equilibrate: integrate until teq
	integrate(r,dr,p,dp,deriv,0,teq,dt);

	for( int n =0;n<navg;n++) {
		integrate(r,dr,p,dp,deriv,0,tf,dt);

	}




	return 0;
}
