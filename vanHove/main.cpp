
#include <iostream>
#include <vector>
#include <string>
#include <random>

#include "../headers/read.h"
#include "../headers/deriv.h"
#include "../headers/integrate.h"
#include "../headers/pair_distribution.h"
#include "../headers/vanHove.h"


using namespace std;


int main(int argc, char *argv[])
{
	// params for vh calc
	int Nt,nbin;
	double tsample;

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
	read_variables(N,rho,Dt,Dr,gamma,beta,eps,sigma,dt,tf,teq,navg,seed,name,input_name);
	L = pow(N/rho,1./3);
	Deriv deriv(N,L,Dt,Dr,gamma,beta,eps,sigma,seed);

	default_random_engine gen(seed);

	vector<vector<double> > pdist_mat(navg);
	vector<vector<double> > r(N,vector<double>(3));
	vector<vector<double> > dr(N,vector<double>(3));
	vector<vector<double> > p(N,vector<double>(3));
	vector<vector<double> > dp(N,vector<double>(3));

	// start with random initial vectors
	rand_vecs(r,N,3,-.5*L,.5*L,gen);
	rand_vecs(p,N,3,-.5*L,.5*L,gen,1.);
	integrate(r,dr,p,dp,deriv,0,teq,dt);

	// calculate the Van Hove function
	// step 1: calculate trajectories
	// save every tsample

	// matrix with al trajectories
	// initialize all trajectories with the current position
	vector<vector<vector<double> > > rt(Nt,r);
	for(int ti = 1;ti<Nt;++ti) {
		integrate(r,dr,p,dp,deriv,0,tsample,dt);
		rt[ti] = r;
	}
	

	// step 2: calculate VH function from trajectories
	vector<vector<double> > vhr(Nt,vector<double>(nbin,0.));
	vanHove_r(rt,vhr);


	return 0;
}




