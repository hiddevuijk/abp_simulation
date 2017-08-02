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
	for(int i=0;i<N;++i) {
		x -= L*round(r[i][xyz]/L);
		j = floor(x/dr);
		pAvg[j][0] += p[i][0]/N;
		pAvg[j][1] += p[i][1]/N;
		pAvg[j][2] += p[i][2]/N;
	}
}



#endif
