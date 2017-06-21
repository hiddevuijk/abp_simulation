
#include <iostream>
#include <vector>
#include <string>

#include "deriv.h"

using namespace std;


int main()
{

	int N = 2;
	double Dt = 1./30;
	double Dr = 1./10.;
	double beta = 1.;
	double eps = 0.;
	double sigma = 0.;
	int seed = 123456789;

	Deriv deriv(N,Dt,Dr,beta,eps,sigma,seed);
	vector<vector<double> > r(N,vector<double>(3,1.));
	vector<vector<double> > dr(N,vector<double>(3));
	vector<vector<double> > p(N,vector<double>(3,0.1));
	vector<vector<double> > dp(N,vector<double>(3));

	for(int i=0;i<2;++i)
		deriv(r,dr,p,dp,0.1);
	cout << r[0][0] << endl;

	return 0;
}





