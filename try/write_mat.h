

#include <fstream>
#include <vector>
#include <string>

void write_mat(const std::vector<std::vector<double> >& m,int nv,
	int nelem, const char outname[], char sep1 = '\n',char sep2 = ';')
{

	std::ofstream out;
	out.open(outname);
	for(int nvi=0;nvi<nv;nvi++){
		for(int nelemi=0;nelemi<nelem;nelemi++){
			out << m[nvi][nelemi];
			if(nelemi + 1 < nelem) out << sep2;
		}
		out << sep1;
	}
}






