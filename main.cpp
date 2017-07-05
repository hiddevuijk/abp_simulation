
#include <iostream>
#include <vector>
#include <string>
#include <random>

#include "deriv.h"

using namespace std;


int main()
{

	int N = 2;
	double Dt = 1./30;
	double Dr = 1./10.;
	double gamma = 1.;
	double beta = 1.;
	double eps = 0.;
	double sigma = 0.;
	double L = 1.;
	int seed = 123456789;

	Deriv deriv(N,L,Dt,Dr,gamma,beta,eps,sigma,seed);

	default_random_engine gen;
	uniform_real_distribution<double> dist(0.,L);
	vector<vector<double> > r(N,vector<double>(3,1.));
	vector<vector<double> > dr(N,vector<double>(3));
	vector<vector<double> > p(N,vector<double>(3,0.1));
	vector<vector<double> > dp(N,vector<double>(3));

	for(int i=0;i<N;i++) {
		r[i][0] = dist(gen); 
		r[i][1] = dist(gen);
		r[i][2] = dist(gen);
		p[i][0] = dist(gen);
		p[i][1] = dist(gen);
		p[i][2] = dist(gen);
	}
	for(int i=0;i<2000;++i){
		deriv(r,dr,p,dp,0.1);
		cout << r[0][0] << endl;
	}
	return 0;
}





