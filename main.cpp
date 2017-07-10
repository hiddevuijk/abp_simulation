
#include <iostream>
#include <vector>
#include <string>
#include <random>

#include "read.h"
#include "deriv.h"
#include "integrate.h"
#include "pair_distribution.h"

using namespace std;


int main(int argc, char *argv[])
{

	int N, seed;
	double Dt,Dr,gamma,beta,eps,sigma,L,dt,tf,teq;
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
	read_variables(N,rho,Dt,Dr,gamma,beta,eps,sigma,dt,tf,teq,seed,name,input_name);
	L = pow(N/rho,1./3);
	Deriv deriv(N,L,Dt,Dr,gamma,beta,eps,sigma,seed);

	default_random_engine gen(seed);

	vector<vector<double> > pdist_mat(nAvg);
	for( int n =0;n<nAvg;n++) {

		vector<vector<double> > r(N,vector<double>(3));
		vector<vector<double> > dr(N,vector<double>(3));
		vector<vector<double> > p(N,vector<double>(3));
		vector<vector<double> > dp(N,vector<double>(3));

		// start with random initial vectors
		rand_vecs(r,N,3,-.5*L,.5*L,gen);
		rand_vecs(p,N,3,-.5*L,.5*L,gen,1.);

		integrate(r,dr,p,dp,deriv,0,teq,dt);

		
		pdist_mat[n] = pair_distances(r,L);
	}
	write_vec(pdist,"pdist.dat");
	return 0;
}




