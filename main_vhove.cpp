
#include <iostream>
#include <vector>
#include <string>
#include <random>

#include "read.h"
#include "deriv.h"
#include "integrate.h"
#include "pair_distribution.h"
#include "diffconst.h"
#include "fcc_lattice.h"
#include "vanHove.h"


using namespace std;


int main(int argc, char *argv[])
{

	// number of particles
	int N;
	//density
	double rho;
	// size of the box
	double L;
	// diff consts
	double Dt,Dr;
	//porential parameters
	double gamma, beta, sigma,eps;
	// Van Hove parameters
	double dt_vh,tm,dr_vh,rm;
	int Nt_vh,Nr_vh;	
	// inital time used to equilibrate
	double teq;
	// integration time step
	double dt;
	//  average over navg vanHove functions
	int navg;
	// seed for the random number generator
	int seed;
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
	read_variables_vhove(N,rho,Dt,Dr,gamma,beta,eps,sigma,dt,dt_vh,tm,dr_vh,rm,teq,navg,seed,name,input_name);
	L = pow(N/rho,1./3);
	Nt_vh = tm/dt_vh+1;
	Nr_vh = rm/dr_vh+1;

	vector<double> t_vh(Nt_vh,0.0);
	vector<double> r_vh(Nr_vh,0.0);
	for(int i=0;i<Nt_vh;++i)
		t_vh[i] = i*dt_vh;
	for(int i=0;i<Nr_vh;++i)
		r_vh[i] = (0.5+i)*dr_vh;

	Deriv deriv(N,L,Dt,Dr,gamma,beta,eps,sigma,seed);

	default_random_engine gen(seed);

	vector<vector<double> > r(N,vector<double>(3));
	vector<vector<double> > dr(N,vector<double>(3));
	vector<vector<double> > p(N,vector<double>(3));
	vector<vector<double> > dp(N,vector<double>(3));

	// start with random p(t=0)
	rand_vecs(p,N,3,-.5*L,.5*L,gen,1.);

	// start with particles on fcc lattic
	// with largest possible lattice spacing
	init_r_fcc(r,N,sigma,L);

	// integrate until teq in order to equilibrate
	integrate(r,dr,p,dp,deriv,0,teq,dt);

	// initial position vectors
	vector<vector<double> > r0 = r;

	// traveld distances at an instance of time, compared to r0
	//vector<double> delta_r(N);

	// Van Hove result
	vector<vector<double> > vhove(Nt_vh,vector<double>(Nr_vh,0.0));

	// van Hove single measurement
	vector<vector<double> > vhove_temp(Nt_vh,vector<double>(Nr_vh,0.0));

	for(int avgi = 0;avgi<navg;++avgi) {
		r0 = r;
		for( int ti =0;ti<Nt_vh;ti++) {
			integrate(r,dr,p,dp,deriv,0,dt_vh,dt);
			get_vhove(vhove_temp[ti],r,r0,dr_vh,Nr_vh);
		}
		add2result(vhove,vhove_temp,navg);
	}
	if(name != "") name += '_';
	write_vec(r_vh,name+"r.dat");
	write_vec(t_vh,name+"t.dat");
	write_mat(vhove,Nt_vh,Nr_vh,name+"vhove.dat");


	return 0;
}




