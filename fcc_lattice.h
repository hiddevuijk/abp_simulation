#ifndef GUARD_fcc_lattice_h
#define GUARD_fcc_lattice_h

#include <vector>
#include <math.h>
#include <iostream> 
void init_r_fcc(std::vector<std::vector<double> >& r,
	int N, double sigma, double L)
{

	// n is the minimum number of basic lattice block needed
	// in order tot host N particles
	int n = ceil(pow(N/3.,1./3));

	// l : lattice const.
	double l = L/n;


	// xi,yi,zi: vector pointing to the point on the lattice block
	// closest to (0,0,0)
	double zi = 0;
	double yi = 0;
	double xi = 0;

	// i: number of placed particles
	int i = 0;
	while(zi<n) {
	 while(yi<n) {
	  while(xi<n) {
		if(i>=N) break;
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
