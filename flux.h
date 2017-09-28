#ifndef GUARD_flux_h
#define GUARD_flux_h

#include <vector>

void flux(
	const std::vector<std::vector<double> >& r1,
	const std::vector<std::vector<double> >& r2,
	std::vector<double>& f,
	int xyz, double dr, double L)
{
	// flux correspond to the right boundary of the bins
	double x1,x2,j1,j2;

	std::fill(f.begin(),f.end(),0.0);

	for(int i=0;i<r1.size();++i) {
		x1 = r1[i][xyz] - L*floor(r1[i][xyz]/L);
		j1 = floor(x1/dr);
		x2 = r2[i][xyz] - L*floor(r2[i][xyz]/L);
		j2 = floor(x2/dr);

		// check |j1-j2|<=1

		if(j2-j1>0)
			f[j1] += 1.;
		if(j2-j1<0)
			f[j2] += 1.;

	}

}




#endif
