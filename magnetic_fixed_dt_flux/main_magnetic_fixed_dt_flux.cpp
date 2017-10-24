
#include <iostream>
#include <vector>
#include <string>
#include <random>

#include "read.h"
#include "deriv_magnetic_fixed_dt.h"
#include "integrate.h"
#include "pair_distribution.h"
#include "fcc_lattice.h"
#include "orientation.h"
#include "density.h"
#include "flux.h"

using namespace std;

typedef vector<double> M1;
typedef vector<vector<double> > M2;
typedef vector<vector<vector<double> > > M3;

int main(int argc, char *argv[])
{

	int navg;
	int N, seed;
	double Dt,Dr,v0,qB,omega0,beta,eps,sigma,L,dt,tf,teq;
	double rho;
	double bs_pAvg;		// binsize for orientation
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
	read_variables_magnetic(N,rho,Dt,Dr,v0,qB,omega0,beta,eps,
			sigma,dt,tf,teq,navg,bs_pAvg,seed,name,input_name);

	// L: box size
	L = pow(N/rho,1./3);

	// omega = omega0 2 pi / L
	double omega = omega0*2*acos(-1)/L;


	int NpAvg = ceil(L/bs_pAvg);
	M2 pAvgx_temp(NpAvg,M1(3,0.0));
	M2 pAvgy_temp(NpAvg,M1(3,0.0));
	M2 pAvgz_temp(NpAvg,M1(3,0.0));
	M2 pAvgx(NpAvg,M1(3,0.0));
	M2 pAvgy(NpAvg,M1(3,0.0));
	M2 pAvgz(NpAvg,M1(3,0.0));
	M1 rhox(NpAvg,0.);
	M1 rhoy(NpAvg,0.);
	M1 rhoz(NpAvg,0.);
	M1 rhox_temp(NpAvg,0.);
	M1 rhoy_temp(NpAvg,0.);
	M1 rhoz_temp(NpAvg,0.);
	M1 rhox_temp1(NpAvg,0.);
	M1 rhoy_temp1(NpAvg,0.);
	M1 rhoz_temp1(NpAvg,0.);
	M1 fluxy(NpAvg,0.);
	M1 fluyx(NpAvg,0.);
	M1 fluyy(NpAvg,0.);
	M1 fluxy_temp(NpAvg,0.);
	M1 fluyx_temp(NpAvg,0.);
	M1 fluyy_temp(NpAvg,0.);
	M3 fxy(NpAvg,M2(NpAvg,M1(2,0.)));


	vector<double> bins(NpAvg,0.0);
	for(int i=0;i<NpAvg;++i)
		bins[i] = (i+.5)*bs_pAvg;


	vector<vector<double> > r(N,vector<double>(3));
	vector<vector<double> > dr(N,vector<double>(3));
	vector<vector<double> > p(N,vector<double>(3,1.));
	vector<vector<double> > dp(N,vector<double>(3,1.));

	// initalize r vector: put particles on a fcc lattice
	init_r_fcc(r,N,sigma,L);

	// init deriv objec to perform integration
	Deriv deriv(N,L,Dt,Dr,v0,qB,omega,beta,eps,sigma,seed);

	// equilibrate: integrate until teq
	integrate(r,dr,p,dp,deriv,0,teq,dt);

	int errorCount = 0;

	for( int n =0;n<navg;n++) {
		integrate(r,dr,p,dp,deriv,0,tf-dt,dt);
		deriv(r,dr,p,dp,dt,0.1,0.1);

		// get the orientation
		orientation(r,p,pAvgx_temp,0,bs_pAvg,L);
		orientation(r,p,pAvgy_temp,1,bs_pAvg,L);
		orientation(r,p,pAvgz_temp,2,bs_pAvg,L);

		rhox_temp = rhox_temp1;
		rhoy_temp = rhoy_temp1;
		rhoz_temp = rhoz_temp1;

		// and the density
		density(r,rhox_temp1,0,bs_pAvg,L);
		density(r,rhoy_temp1,1,bs_pAvg,L);
		density(r,rhoz_temp1,2,bs_pAvg,L);

		fluxXY(r,dr,fxy,L,bs_pAvg);

		flux2(r,dr,fluxy_temp,0,1,bs_pAvg,L,errorCount);
		flux2(r,dr,fluyx_temp,1,0,bs_pAvg,L,errorCount);
		flux2(r,dr,fluyy_temp,1,1,bs_pAvg,L,errorCount);


		for(int i=0;i<NpAvg;++i) {
			pAvgx[i][0] += pAvgx_temp[i][0]/(navg*N*bs_pAvg);
			pAvgx[i][1] += pAvgx_temp[i][1]/(navg*N*bs_pAvg);
			pAvgx[i][2] += pAvgx_temp[i][2]/(navg*N*bs_pAvg);
			
			pAvgy[i][0] += pAvgy_temp[i][0]/(navg*N*bs_pAvg);
			pAvgy[i][1] += pAvgy_temp[i][1]/(navg*N*bs_pAvg);
			pAvgy[i][2] += pAvgy_temp[i][2]/(navg*N*bs_pAvg);
			
			pAvgz[i][0] += pAvgz_temp[i][0]/(navg*N*bs_pAvg);
			pAvgz[i][1] += pAvgz_temp[i][1]/(navg*N*bs_pAvg);
			pAvgz[i][2] += pAvgz_temp[i][2]/(navg*N*bs_pAvg);

			rhox[i] += rhox_temp[i]/(navg*N*bs_pAvg);
			rhoy[i] += rhoy_temp[i]/(navg*N*bs_pAvg);
			rhoz[i] += rhoz_temp[i]/(navg*N*bs_pAvg);

			fluxy[i] += fluxy_temp[i]/(navg*dt*L*L);
			fluyx[i] += fluyx_temp[i]/(navg*dt*L*L);
			fluyy[i] += fluyy_temp[i]/(navg*dt*L*L);


		}	
	}

	for(int i=0;i<NpAvg;++i) {
		for(int j=0;j<NpAvg;++j) {
			fxy[i][j][0] /= navg*N*dt*bs_pAvg*bs_pAvg;
			fxy[i][j][1] /= navg*N*dt*bs_pAvg*bs_pAvg;
		}
	}



	
	
	write_mat(pAvgx,pAvgx.size(),pAvgx[0].size(),"pAvgx.dat");
	write_mat(pAvgy,pAvgy.size(),pAvgy[0].size(),"pAvgy.dat");
	write_mat(pAvgz,pAvgz.size(),pAvgz[0].size(),"pAvgz.dat");
	write_vec(rhox,"rhox.dat");
	write_vec(rhoy,"rhoy.dat");
	write_vec(rhoz,"rhoz.dat");
	write_vec(fluxy,"fluxy.dat");
	write_vec(fluyx,"fluyx.dat");
	write_vec(fluyy,"fluyy,dat");
	write_vec(bins,"bins.dat");
	

	ofstream outX;
	outX.open("fxyX.dat");
	ofstream outY;
	outY.open("fxyY.dat");

	for(int n=0;n<fxy.size();n++){
		for(int m=0;m<fxy[n].size();m++){
			outX << fxy[n][m][0];
			outY << fxy[n][m][1];
			if(m + 1 < fxy[n].size()) {
				outX << ';';
				outY << ';';
			}
		}
		outX << '\n';
		outY << '\n';
	}




	cout << errorCount;

	return 0;
}
