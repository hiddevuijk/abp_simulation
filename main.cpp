
#include <iostream>
#include <vector>
#include <string>
#include <random>

#include "deriv.h"
#include "integrate.h"
using namespace std;


int main()
{

	int N = 100;
	double Dt = 1./1;
	double Dr = 1./10.;
	double gamma = 1.;
	double beta = 1.;
	double eps = 1.;
	double sigma = 1.;
	double L = 10.;
	int seed = 123456789;
	double dt = (1.e-5)*sigma*sigma/Dt;
	double tf = 1e1;
	double rho = 0.001;
	L = pow(N/rho,1./3);
	cout << "L  = " << L << endl;
	cout << "dt = " << dt << endl << endl;
	

	Deriv deriv(N,L,Dt,Dr,gamma,beta,eps,sigma,seed);

	default_random_engine gen(seed);
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
	cout << r[0][0] << endl;

	integrate(r,dr,p,dp,deriv,0,tf,dt);
	cout << r[0][0] << endl;

	return 0;
}





