#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

#include "write_mat.h"
using namespace std;

int main()
{
	int nv = 10;
	int lv = 3;

	vector<vector<double> > M(nv,vector<double>(3,0));
	
	vector<double> vec(lv);
	for(int nvi=0;nvi<nv;nvi++){
		for(int lvi=0;lvi<lv;lvi++) {
			M[nvi][lvi] = nvi + lvi/100.;
		}
	}

	write_mat(M,nv,lv,"m.dat");	




	return 0;
}



