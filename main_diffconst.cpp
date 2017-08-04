
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

using namespace std;


int main(int argc, char *argv[])
{

	int Nt;
	int N, seed;
	double Dt,Dr,v0,gamma,beta,eps,sigma,L,dt,tf,teq;
	double rho;
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
	read_variables_diffconst(N,rho,Dt,Dr,v0,
		gamma,beta,eps,sigma,
		dt,tf,teq,Nt,seed,name,input_name);
	L = pow(N/rho,1./3);

	// -1 -> v0 = v
	Deriv deriv(N,L,Dt,Dr,v0,-1,gamma,beta,eps,sigma,seed);

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
	vector<double> delta_r(N);

	// <(r-r0)^2>, avg over all particles
	vector<double> dR(Nt);

	for( int ti =0;ti<Nt;ti++) {
		integrate(r,dr,p,dp,deriv,0,tf,dt);

		// calculate <(r(t=ti*tf)-r(t=t0))^2> 
		dR[ti] = get_dR(r,r0);
	}

	write_vec(dR,name+".dat",Nt);	


	return 0;
}




