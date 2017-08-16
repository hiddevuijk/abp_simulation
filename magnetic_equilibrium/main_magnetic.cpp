
#include <iostream>
#include <vector>
#include <string>
#include <random>

#include "read.h"
#include "deriv_magnetic.h"
#include "integrate.h"
#include "pair_distribution.h"
#include "fcc_lattice.h"
#include "orientation.h"
#include "density.h"

using namespace std;


int main(int argc, char *argv[])
{

	int navg;
	int N, seed;
	double Dt,qB,omega0,gamma,L,dt,tf,teq;
	double rho;
	double bs;
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
	read_variables_magnetic(N,rho,Dt,qB,omega0,gamma,
			dt,tf,teq,navg,bs,seed,name,input_name);
	// L: box size
	L = pow(N/rho,1./3);

	// omega = omega0 2 pi / L
	double omega = omega0*2*acos(-1)/L;

	int Nbin = ceil(L/bs);
	vector<double> rhox(Nbin,0.);
	vector<double> rhoy(Nbin,0.);
	vector<double> rhoz(Nbin,0.);
	vector<double> rhox_temp(Nbin,0.);
	vector<double> rhoy_temp(Nbin,0.);
	vector<double> rhoz_temp(Nbin,0.);


	vector<double> bins(Nbin,0.0);
	for(int i=0;i<Nbin;++i)
		bins[i] = (i+0.5)*bs;


	vector<vector<double> > r(N,vector<double>(3));
	vector<vector<double> > dr(N,vector<double>(3));

	// initalize r vector: put particles on a fcc lattice
	init_r_fcc(r,N,L);

	// init deriv objec to perform integration
	Deriv deriv(N,L,Dt,qB,omega,gamma,seed);

	// equilibrate: integrate until teq
	integrate_nof(r,dr,deriv,0,teq,dt);
	for( int n =0;n<navg;n++) {
		integrate(r,dr,deriv,0,tf,dt);

		density(r,rhox_temp,0,bs,L);
		density(r,rhoy_temp,1,bs,L);
		density(r,rhoz_temp,2,bs,L);

		for(int i=0;i<Nbin;++i) {
			rhox[i] += rhox_temp[i]/(navg*N);
			rhoy[i] += rhoy_temp[i]/(navg*N);
			rhoz[i] += rhoz_temp[i]/(navg*N);
		}	
	}

	
	write_vec(rhox,"rhox.dat");
	write_vec(rhoy,"rhoy.dat");
	write_vec(rhoz,"rhoz.dat");
	write_vec(bins,"bins.dat");



	return 0;
}
