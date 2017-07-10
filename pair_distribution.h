#ifndef GUARD_pair_distribution_h
#define GUARD_pair_distribution_h

#include <vector>
#include <algorithm>
#include <math.h>
#include "vecmanip.h"

namespace pair_distribution_vars {
	const double PI4 = 4*std::acos(-1.);
}

std::vector<double> pair_distribution(
	const std::vector<std::vector<double> >& r,
	double L, int nbin)
{
	int N = r.size();
	std::vector<double> dist(N*(N-1)/2,0.);
	int di =0;
	for(int i=0;i<N;++i) {
		for(int j=i+1;j<N;++j) {
			dist[di++] = distance_3d(r[i],r[j],L);
		}
	}	

	std::sort(dist.begin(),dist.end());
	
	std::vector<double> pd_vec(nbin);
	di =0;
	double bs = L/nbin;
	for(int i=0;i<nbin;++i)
		while(dist[di++]<(i+1)*bs)
			pd_vec[i] += 1.;

	// normalize with density and volume of the shell
	for(int i=0;i<nbin;++i) {
		pd_vec[i] /= N;
		pd_vec[i] *=(L*L*L);
		pd_vec[i] /= pair_distribution_vars::PI4*(i+0.5)*bs*(i+0.5)*bs*bs;
	}
	return pd_vec;
}



#endif 
