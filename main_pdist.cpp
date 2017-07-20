
#include <iostream>
#include <vector>
#include <string>
#include <random>

#include "read.h"
#include "deriv.h"
#include "integrate.h"
#include "pair_distribution.h"
#include "fcc_lattice.h"


using namespace std;


int main(int argc, char *argv[])
{

	int navg;
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
	read_variables_pairdist(N,rho,Dt,Dr,gamma,beta,eps,sigma,dt,tf,teq,navg,seed,name,input_name);
	L = pow(N/rho,1./3);
	int ndist = navg*N*(N-1)/2;

	vector<double> pdist(ndist);
	vector<double> pdist_temp(N*(N-1)/2);
	vector<vector<double> > r(N,vector<double>(3));
	vector<vector<double> > dr(N,vector<double>(3));
	vector<vector<double> > p(N,vector<double>(3));
	vector<vector<double> > dp(N,vector<double>(3));

	// initalize r vector: put particles on a fcc lattice
	init_r_fcc(r,N,sigma,L);

	// init deriv objec to perform integration
	Deriv deriv(N,L,Dt,Dr,gamma,beta,eps,sigma,seed);

	// equilibrate: integrate until teq
	integrate(r,dr,p,dp,deriv,0,teq,dt);

	for( int n =0;n<navg;n++) {
		integrate(r,dr,p,dp,deriv,0,tf,dt);
		pdist_temp = pair_distances(r,L);
		for(int i =0;i<N*(N-1)/2;++i)
			pdist[i+n*N*(N-1)/2] = pdist_temp[i];
	}
	std::sort(pdist.begin(),pdist.end());
	write_vec(pdist,"pdist.dat");
	return 0;
}
