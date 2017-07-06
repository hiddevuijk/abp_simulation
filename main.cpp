
#include <iostream>
#include <vector>
#include <string>
#include <random>

#include "deriv.h"
#include "integrate.h"
#include "read.h"


using namespace std;


int main()
{

	int N, seed;
	double Dt,Dr,gamma,beta,eps,sigma,L,dt,tf,teq;

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
	read_variables(N,Dt,Dr,gamma,beta,eps,sigma,L,dt,tf,teq,seed,name,input_name);

	Deriv deriv(N,L,Dt,Dr,gamma,beta,eps,sigma,seed);

	default_random_engine gen(seed);

	vector<vector<double> > r(N,vector<double>(3));
	vector<vector<double> > dr(N,vector<double>(3));
	vector<vector<double> > p(N,vector<double>(3));
	vector<vector<double> > dp(N,vector<double>(3));

	// start with random initial vectors
	rand_vecs(r,N,3,-.5*L,.5*L,gen);
	rand_vecs(p,N,3,-.5*,.5*L,gen,1.);

	integrate(r,dr,p,dp,deriv,0,tf,dt);
	cout << r[0][0] << endl;

	return 0;
}





