
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
		t_vh[i] = (i+1)*dt_vh;
	// r_vh contains the right values of the bins
	for(int i=0;i<Nr_vh;++i)
		r_vh[i] = (i+1)*dr_vh;

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


	// Van Hove in the three directions result
	vector<vector<double> > vhovex(Nt_vh,vector<double>(Nr_vh,0.0));
	vector<vector<double> > vhovey(Nt_vh,vector<double>(Nr_vh,0.0));
	vector<vector<double> > vhovez(Nt_vh,vector<double>(Nr_vh,0.0));
	
	// van Hove single measurement
	vector<vector<double> > vhovex_temp(Nt_vh,vector<double>(Nr_vh,0.0));
	vector<vector<double> > vhovey_temp(Nt_vh,vector<double>(Nr_vh,0.0));
	vector<vector<double> > vhovez_temp(Nt_vh,vector<double>(Nr_vh,0.0));


	for(int avgi = 0;avgi<navg;++avgi) {
		r0 = r;
		for( int ti =0;ti<Nt_vh;ti++) {
			integrate(r,dr,p,dp,deriv,0,dt_vh,dt);
			get_vhovexyz(vhovex_temp[ti],vhovey_temp[ti],vhovez_temp[ti],r,r0,dr_vh,Nr_vh);
		}
		add2result(vhovex,vhovex_temp,navg);
		add2result(vhovey,vhovey_temp,navg);
		add2result(vhovez,vhovez_temp,navg);
	}

	// shift r_vh with -0.5*dr_vh to get mid points of the bins
	for(int i=0;i<r_vh.size();++i)
		r_vh[i] -= 0.5*dr_vh;


	if(name != "") name += '_';
	write_vec(r_vh,name+"r.dat");
	write_vec(t_vh,name+"t.dat");
	write_mat(vhovex,Nt_vh,Nr_vh,name+"vhove_x.dat");
	write_mat(vhovey,Nt_vh,Nr_vh,name+"vhove_y.dat");
	write_mat(vhovez,Nt_vh,Nr_vh,name+"vhove_z.dat");



	return 0;
}




