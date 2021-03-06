#ifndef GUARD_flux_h
#define GUARD_flux_h

#include <vector>
#include <assert.h>
#include <iostream>
#include <math.h>

namespace Flux {
	typedef std::vector<std::vector<double> > M2;
	typedef std::vector<std::vector<std::vector<double> > > M3;
};

void fluxXY( const Flux::M2& r,const Flux::M2& v,
	Flux::M3& f,double L, double bs)
{
	double x,y;
	int xi,yi;
	for(int i=0;i<r.size();++i) {
		x = r[i][0] - L*floor(r[i][0]/L);
		y = r[i][1] - L*floor(r[i][1]/L);
		xi = floor(x/bs);
		yi = floor(y/bs);
		f[xi][yi][0] += v[i][0];
		f[xi][yi][1] += v[i][1];
	} 

}


void flux2(
	const std::vector<std::vector<double> >& r,
	const std::vector<std::vector<double> >& dr,
	std::vector<double>& f,
	int fxyz, int xyz, double bs, double L,int errorCount)
{
	// flux correspond to the right boundary of the bins
	double x1,x2,r1ixyz;
	int j1,j2,dj;
	std::fill(f.begin(),f.end(),0.0);

	for(int i=0;i<r.size();++i) {
		x2 = r[i][fxyz] - L*floor(r[i][fxyz]/L);
		j2 = floor(x2/bs);
		r1ixyz = r[i][fxyz] - dr[i][fxyz];
		x1 = r1ixyz - L*floor(r1ixyz/L);

		j1 = floor(x1/bs);

		dj = fabs(j2-j1);
		dj = (dj)%(f.size()-1);
		if(dj>1) {
			dj = 1;
			++errorCount;
		}
		if(j2<j1) dj *= -1;
	
		//assert(dj<=1 and dj>=-1);
		r1ixyz = r[i][xyz] - dr[i][xyz];
		x1 = r1ixyz - L*floor(r1ixyz/L);
		j1 = floor(x1/bs);
		if(dj==1)
			f[j1] += 1.;
		if(dj==-1)
			f[j1] -= 1.;

	}

}



void flux(
	const std::vector<std::vector<double> >& r,
	const std::vector<std::vector<double> >& dr,
	std::vector<double>& f,
	int xyz, double bs, double L)
{
	// flux correspond to the right boundary of the bins
	double x1,x2,r1ixyz;
	int j1,j2,dj;
	std::fill(f.begin(),f.end(),0.0);

	for(int i=0;i<r.size();++i) {
		x2 = r[i][xyz] - L*floor(r[i][xyz]/L);
		j2 = floor(x2/bs);
		r1ixyz = r[i][xyz] - dr[i][xyz];
		x1 = r1ixyz - L*floor(r1ixyz/L);
		j1 = floor(x1/bs);
		dj = fabs(j2-j1);
		dj = (dj)%(f.size()-1);
		if(j2<j1) dj *= -1;
		// check |j1-j2|<=1
		assert(dj<=1 and dj>=-1);
		if(dj==1)
			f[j1] += 1.;
		if(dj==-1)
			f[j2] -= 1.;

	}

}




#endif
