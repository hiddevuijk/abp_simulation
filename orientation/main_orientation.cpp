
#include <iostream>
#include <vector>
#include <string>
#include <random>

#include "read.h"
#include "deriv.h"
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
	double Dt,Dr,v0,omega0,gamma,beta,eps,sigma,L,dt,tf,teq;
	double rho;
	double bs_pAvg;		// binsize for orientation
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
	read_variables_orientation(N,rho,Dt,Dr,v0,omega0,gamma,beta,eps,sigma,dt,tf,teq,navg,bs_pAvg,seed,name,input_name);


	L = pow(N/rho,1./3);

	// omega = omega0 2 pi / L
	double omega = omega0*2*acos(-1)/L;

	int NpAvg = ceil(L/bs_pAvg);
	vector<vector<double> > pAvgx_temp(NpAvg,vector<double>(3,0.0));
	vector<vector<double> > pAvgy_temp(NpAvg,vector<double>(3,0.0));
	vector<vector<double> > pAvgz_temp(NpAvg,vector<double>(3,0.0));
	vector<vector<double> > pAvgx(NpAvg,vector<double>(3,0.0));
	vector<vector<double> > pAvgy(NpAvg,vector<double>(3,0.0));
	vector<vector<double> > pAvgz(NpAvg,vector<double>(3,0.0));
	vector<double> rhox(NpAvg,0.);
	vector<double> rhoy(NpAvg,0.);
	vector<double> rhoz(NpAvg,0.);
	vector<double> rhox_temp(NpAvg,0.);
	vector<double> rhoy_temp(NpAvg,0.);
	vector<double> rhoz_temp(NpAvg,0.);


	vector<double> pAvg_bins(NpAvg,0.0);

	for(int i=0;i<NpAvg;++i)
		pAvg_bins[i] = (i+1)*bs_pAvg;


	vector<vector<double> > r(N,vector<double>(3));
	vector<vector<double> > dr(N,vector<double>(3));
	vector<vector<double> > dp(N,vector<double>(3,1.));
	vector<vector<double> > p(N,vector<double>(3,0.));

	//random initial orientation
	default_random_engine generator(seed-1);
	normal_distribution<double> ndist(0.,1.);
	for(int i=0;i<N;++i) {
		p[i][0] = ndist(generator);
		p[i][1] = ndist(generator);
		p[i][2] = ndist(generator);
		normalize(p[i]);
	}

	// initalize r vector: put particles on a fcc lattice
	init_r_fcc(r,N,sigma,L);

	// init deriv objec to perform integration
	Deriv deriv(N,L,Dt,Dr,v0,omega,gamma,beta,eps,sigma,seed);

	// equilibrate: integrate until teq
	integrate(r,dr,p,dp,deriv,0,teq,dt);

	for( int n =0;n<navg;n++) {
		integrate(r,dr,p,dp,deriv,0,tf,dt);

		// get the orientation
		orientation(r,p,pAvgx_temp,0,bs_pAvg,L);
		orientation(r,p,pAvgy_temp,1,bs_pAvg,L);
		orientation(r,p,pAvgz_temp,2,bs_pAvg,L);

		density(r,rhox_temp,0,bs_pAvg,L);
		density(r,rhoy_temp,1,bs_pAvg,L);
		density(r,rhoz_temp,2,bs_pAvg,L);


		for(int i=0;i<NpAvg;++i) {
			pAvgx[i][0] += pAvgx_temp[i][0]/navg;
			pAvgx[i][1] += pAvgx_temp[i][1]/navg;
			pAvgx[i][2] += pAvgx_temp[i][2]/navg;
			
			pAvgy[i][0] += pAvgy_temp[i][0]/navg;
			pAvgy[i][1] += pAvgy_temp[i][1]/navg;
			pAvgy[i][2] += pAvgy_temp[i][2]/navg;
			
			pAvgz[i][0] += pAvgz_temp[i][0]/navg;
			pAvgz[i][1] += pAvgz_temp[i][1]/navg;
			pAvgz[i][2] += pAvgz_temp[i][2]/navg;

			rhox[i] += rhox_temp[i]/navg;
			rhoy[i] += rhoy_temp[i]/navg;
			rhoz[i] += rhoz_temp[i]/navg;


		}	
	}


	// shift bins with .5bs to mid point
	for(int i=0;i<NpAvg;++i)
		pAvg_bins[i] -= 0.5*bs_pAvg;
	
	write_mat(pAvgx,pAvgx.size(),pAvgx[0].size(),"pAvgx.dat");
	write_mat(pAvgy,pAvgy.size(),pAvgy[0].size(),"pAvgy.dat");
	write_mat(pAvgz,pAvgz.size(),pAvgz[0].size(),"pAvgz.dat");
	write_vec(pAvg_bins,"pAvg_bins.dat");

	write_vec(rhox,"rhox.dat");
	write_vec(rhoy,"rhoy.dat");
	write_vec(rhoz,"rhoz.dat");


	return 0;
}
