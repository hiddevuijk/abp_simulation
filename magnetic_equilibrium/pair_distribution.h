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
	int j;
	for(int i=0;i<Npd;++i) {
		j = floor(pd[i]/bs);
		if(j<Ng) g[j] += (N-1.)/Npd;
		else g[Ng-1] += (N-1.)/Npd;
	}

	return g;
}


void normalize_g(std::vector<double>& g,
	const std::vector<double>& gr,
	double rho, int N)
{
	for(int i=0;i<g.size();++i){
		// devide by rho*4*pi*r^2*dr
		g[i] /= rho*4*PD::pi*gr[i]*gr[i]*(gr[1]-gr[0]);
	}
}





#endif 
