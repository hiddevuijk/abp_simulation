#ifndef GUARD_pair_distribution_h
#define GUARD_pair_distribution_h

#include <vector>
#include <algorithm>
#include <math.h>
#include "vecmanip.h"

#include <iostream>

namespace PD {
	double pi = acos(-1);
}


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
	std::sort(dist.begin(),dist.end());
	
	return dist;
}


std::vector<double> pd2g(
	const std::vector<double>& pd,
	const std::vector<double>& gr,
	int N,double rho)
{
	int Npd = pd.size();
	int Ng = gr.size();
	std::vector<double> g(Ng,0.0);
	double bs = gr[1] - gr[0];
	int pdi = 0;
	for(int gi=0;gi<Ng;++gi) {
		while(pd[pdi]< gr[gi] ) {
			g[gi] += 1;
			pdi += 1;
			if(pdi>=Npd) break;
		}
		if(pdi>=Npd) break;
	}
	for(int i=0;i<Ng;++i) {
		// devide by 4 pi r^2
		g[i] /= rho*4*PD::pi*(gr[i]-bs)*(gr[i]-bs)*bs;
		g[i] /= Npd;
		g[i] *= N;
	}
	return g;
}








#endif 
