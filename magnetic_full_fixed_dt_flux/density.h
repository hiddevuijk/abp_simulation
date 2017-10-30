#ifndef GUARD_density_h
#define GUARD_density_h


#include <vector>
#include <iostream>
void density(const std::vector<std::vector<double> >& r,
	std::vector<double>& rhoAvg,
	int xyz, double dr,double L)
{
	int N = r.size();
	int j;
	double x;


	// set rhoAvg to zero
	std::fill(rhoAvg.begin(),rhoAvg.end(),0.);


	for(int i=0;i<N;++i) {
		x = r[i][xyz] -  L*floor(r[i][xyz]/L);
		j = floor(x/dr);
		rhoAvg[j] += 1.;
	}
}



#endif
