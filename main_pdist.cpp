
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
	double Dt,Dr,v0,gamma,beta,eps,sigma,L,dt,tf,teq;
	double rho;
	double bs;	// binsize for g(r)
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
	read_variables_pairdist(N,rho,Dt,Dr,v0,gamma,beta,eps,sigma,dt,tf,teq,navg,bs,seed,name,input_name);
	L = pow(N/rho,1./3);


	// vec with distances between particles at snapshot time
	vector<double> pdist_temp(N*(N-1)/2);


	int Ngr = ceil(L/bs);
	vector<double> g(Ngr,0.0);
	vector<double> g_temp(Ngr,0.0);
	// gr contains the 'x' values of g(r), shift -bs in the end
	// to get centers of bins.
	vector<double> gr(Ngr);
	for(int i=0;i<Ngr;++i)
		gr[i] = (i+1)*bs;

	vector<vector<double> > r(N,vector<double>(3));
	vector<vector<double> > dr(N,vector<double>(3));
	vector<vector<double> > p(N,vector<double>(3,1.));
	vector<vector<double> > dp(N,vector<double>(3,1.));

	// initalize r vector: put particles on a fcc lattice
	init_r_fcc(r,N,sigma,L);

	// init deriv objec to perform integration
	// omega = -1 -> v0 = v
	Deriv deriv(N,L,Dt,Dr,v0,-1,gamma,beta,eps,sigma,seed);

	// equilibrate: integrate until teq
	integrate(r,dr,p,dp,deriv,0,teq,dt);



	for( int n =0;n<navg;n++) {
		integrate(r,dr,p,dp,deriv,0,tf,dt);
		pdist_temp = pair_distances(r,L);
		g_temp = pd2g(pdist_temp,gr,N,rho);
		for(int i=0;i<Ngr;++i)
			g[i] += g_temp[i]/navg;	
	}

	// shift gr with -bs to get centers of bins
	for(int i=0;i<Ngr;++i)
		gr[i] -= 0.5*bs;
	normalize_g(g,gr,rho,N);
	write_vec(g,"g.dat");
	write_vec(gr,"gr.dat");

	return 0;
}
