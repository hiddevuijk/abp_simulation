#ifndef GUARD_orientation_h
#define GUARD_orientation_h


#include <vector>

void orientation(const std::vector<std::vector<double> >& r,
	const std::vector<std::vector<double> >& p,
	std::vector<std::vector<double> >& pAvg,
	int xyz, double dr,double L)
{
	int N = r.size();
	int j;
	double x;

	// keeps track of the number of particles in a bin
	std::vector<double> pAvgN(pAvg.size(),0.0);

	// set pAvg to zero
	for(int i=0;i<pAvg.size();++i) {
		pAvg[i][0] = 0;
		pAvg[i][1] = 0;
		pAvg[i][2] = 0;
	}


	for(int i=0;i<N;++i) {

		x = r[i][xyz] -  L*floor(r[i][xyz]/L);
		j = floor(x/dr);
		pAvg[j][0] += p[i][0];
		pAvg[j][1] += p[i][1];
		pAvg[j][2] += p[i][2];

		++pAvgN[j];
	}
	//Divide pAvg[i] by the number of particles in the bin
	for(int i=0;i<pAvg.size();++i) {
		pAvg[i][0] /= pAvgN[i];
		pAvg[i][1] /= pAvgN[i];
		pAvg[i][2] /= pAvgN[i];
	}
}



#endif
