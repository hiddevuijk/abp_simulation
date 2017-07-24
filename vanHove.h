#ifndef GUARD_vanHove_h
#define GUARD_vanHove_h

#include <vector>
#include <math.h>



namespace VH{
	double pi = acos(-1.);
}

void get_vhove(std::vector<double>& vhove_ti,
	const std::vector<std::vector<double> >& r,
	const std::vector<std::vector<double> >& r0,
	double dr_vh, int Nr_vh)
{

	double d;
	int N = r.size();
	int j;
	for(int i=0;i<N;++i) {
		d = distance_3d(r[i],r0[i]);
		j = floor(d/dr_vh);
		if(j<Nr_vh) vhove_ti[j] += 1./N;
		else vhove_ti[Nr_vh-1] += 1./N;

	}



}

void add2result(std::vector<std::vector<double> >& vhove,
	const std::vector<std::vector<double> >& temp,
	int navg)
{

	int Nr = vhove.size();
	int Nt = vhove[0].size();
	for(int i=0;i<Nr;++i) {
		for(int j=0;j<Nt;++j) {
			vhove[i][j] += temp[i][j]/navg;
		}
	}
}

void normalize_vh(std::vector<std::vector<double> >& vhove,
	const std::vector<double>& r)
{
	double dr = r[1] - r[0];
	int Nt_vh = vhove.size();
	int Nr_vh = vhove[0].size();
	double norm;
	for(int ti = 0;ti<Nt_vh;++ti) {
		for(int ri=0;ri<Nr_vh;++ri) {
			norm = r[ri]*r[ri]*dr*4*VH::pi;
			vhove[ti][ri] /= norm;
		}
	}
}









#endif
