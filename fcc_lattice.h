#ifndef GUARD_fcc_lattice_h
#define GUARD_fcc_lattice_h

#include <vector>
#include <math.h>

void init_r_fcc(std::vector<std::vector<double> >& r,
	int N, double sigma, double L)
{

	double l = sqrt(2)*sigma;
	int n = floor(L/l);
	// if (n-1)^3 < N ???
	double zi = 0;
	double yi = 0;
	double xi = 0;
	int i = 0;
	while(zi<n) {
	 while(yi<n) {
	  while(xi<n) {
	  	r[i][0] = xi*l;
		r[i][1] = (yi+0.5)*l;
		r[i++][2] = (zi+0.5)*l;
		if(i>=N) break;
		
	  	r[i][0] = (xi+0.5)*l;
		r[i][1] = yi*l;
		r[i++][2] = (zi+0.5)*l;
		if(i>=N) break;
		
	  	r[i][0] = (xi+0.5)*l;
		r[i][1] = (yi+0.5)*l;
		r[i++][2] = zi*l;
		if(i>=N) break;
		++xi;
	  }

	  if(i>=N) break;
	  xi = 0;
	  ++yi;
	 }

	 if(i>=N) break;
	 yi = 0;
	 ++zi;
	}
}


#endif
