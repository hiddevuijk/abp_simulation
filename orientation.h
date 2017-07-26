#ifndef GUARD_orientation_h
#define GUARD_orientation_h


#include <vector>


std::vector<double> orientation(const std::vector<std::vector<double> >& p)
{

	int N = p.size();
	std::vector<double> pAvg(3,0.0);
	for(int i=0;i<N;++i) {
		pAvg[0] += p[i][0];
		pAvg[1] += p[i][1];
		pAvg[2] += p[i][2];
	}
	pAvg[0] /= N;
	pAvg[1] /= N;
	pAvg[2] /= N;

	return pAvg;
}








#endif
