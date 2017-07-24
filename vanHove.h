#ifndef GUARD_vanHove_h
#define GUARD_vanHove_h

#include <vector>
#include <math.h>


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



#endif
