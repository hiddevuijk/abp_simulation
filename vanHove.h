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

	// clear vhove_ti
	std::fill(vhove_ti.begin(),vhove_ti.end(),0.0);

	for(int i=0;i<N;++i) {
		d = distance_3d(r[i],r0[i]);
		j = floor(d/dr_vh);
		if(j<Nr_vh) vhove_ti[j] += 1./N;
		else vhove_ti[Nr_vh-1] += 1./N;

	}



}

void get_vhovexyz(std::vector<double>& vhovex_ti,
	std::vector<double>& vhovey_ti,
	std::vector<double>& vhovez_ti,
	const std::vector<std::vector<double> >& r,
	const std::vector<std::vector<double> >& r0,
	double dr_vh, int Nr_vh)
{

	double dx,dy,dz;
	int N = r.size();
	int jx,jy,jz;

	// clear vhove{x,y,z}_ti
	std::fill(vhovex_ti.begin(),vhovex_ti.end(),0.0);
	std::fill(vhovey_ti.begin(),vhovey_ti.end(),0.0);
	std::fill(vhovez_ti.begin(),vhovez_ti.end(),0.0);

	for(int i=0;i<N;++i) {
		dx = fabs(r[i][0]-r0[i][0]);
		jx = floor(dx/dr_vh);
		if(jx<Nr_vh) vhovex_ti[jx] += 1./N;
		else vhovex_ti[Nr_vh-1] += 1./N;

		dy = fabs(r[i][1]-r0[i][1]);
		jy = floor(dy/dr_vh);
		if(jy<Nr_vh) vhovey_ti[jy] += 1./N;
		else vhovey_ti[Nr_vh-1] += 1./N;

		dz = fabs(r[i][2]-r0[i][2]);
		jz = floor(dz/dr_vh);
		if(jz<Nr_vh) vhovez_ti[jz] += 1./N;
		else vhovez_ti[Nr_vh-1] += 1./N;

	}

}
void add2result(std::vector<std::vector<double> >& vhove,
	const std::vector<std::vector<double> >& temp,
	int navg)
{

	int Nt = vhove.size();
	int Nr = vhove[0].size();
	for(int i=0;i<Nt;++i) {
		for(int j=0;j<Nr;++j) {
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
