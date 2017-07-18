#ifndef GUARD_pair_distribution_h
#define GUARD_pair_distribution_h

#include <vector>
#include <algorithm>
#include <math.h>
#include "vecmanip.h"

#include <iostream>

namespace pair_dist_vars {
	const double PI4 = 4*std::acos(-1.);
}

std::vector<double> pair_distances(
	const std::vector<std::vector<double> >& r,double L)
{

	int N = r.size();
	std::vector<double> dist(N*(N-1)/2,0.);
	int di =0;
	
	for(int i=0;i<N;++i) {
		for(int j=i+1;j<N;++j) {
			dist[di++] = distance_3d(r[i],r[j],L);
		}
	}	
	//std::sort(dist.begin(),dist.end());
	
	return dist;
}



#endif 
