#ifndef GUARD_vecmanip_h
#define GUARD_vecmanip_h

/* 
 * some functions to manipulate std vectors
*/


#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <random>


double distance_3d(const std::vector<double>& a,
	const std::vector<double>& b, double L=-1)
{
	double d = 0;
	double di = 0;
	if( L < 0 ) {
		for(int i=0;i<3;++i)
			d+= (a[i]-b[i])*(a[i]-b[i]);
	} else {
		for(int i=0;i<3;++i) {
			di = a[i] - b[i];
			di -= L*round(di/L);
			d += di*di;
		}
	}
	return sqrt(d);
}

double dist_sq(const std::vector<double>& a,
	const std::vector<double>& b)
{
	double d = 0;
	for(int i=0;i<3;++i)
		d += (a[i]-b[i])*(a[i]-b[i]);
	return d;
}

// multiply all elements of a vec with a const
void mult_by(std::vector<double>& v, double a)
{
	for(int i=0;i<v.size();++i)
		v[i] *= a;
}

std::vector<double> mult(const std::vector<double>& v, double a)
{
	std::vector<double> out = v;
	mult_by(out,a);
	return out;
}

// returns the length of a vector
double len_vec(const std::vector<double>& a)
{
	double s =0;
	int n = a.size();
	for( int i=0; i<n;++i)
		s += a[i]*a[i];
	return std::sqrt(s);
}

// normalize a vector to length m
// default m=1/
void normalize(std::vector<double>& a,double m = 1.)
{
	double l = len_vec(a);
	int n = a.size();
	for(int i=0;i<n;++i) {
		a[i] /= l;
		a[i] *= m;
	}
}

// create a new vector with entries ci = ai + bi
std::vector<double> add(const std::vector<double>& a,
		const std::vector<double>& b)
{
	int s = a.size();
	std::vector<double> c(s);
	for( int i=0;i<s;++i)
		c[i] = a[i] + b[i];
	return c;

}

// add vector b to vector a
void add_to(std::vector<double>& a,
		const std::vector<double>& b)
{
	int s = a.size();
	for( int i=0;i<s;++i)
		a[i] += b[i];

}

// create new vector with entries ci = ai - bi
std::vector<double> subtr(const std::vector<double>& a,
		const std::vector<double>& b)
{
	int s = a.size();
	std::vector<double> c(a);
	for( int i=0;i<s;++i)
		c[i] = a[i] - b[i];
	return c;

}

// subtract vector b from a
void subtr_of(std::vector<double>& a,
		const std::vector<double>& b)
{
	int s = a.size();
	for( int i=0;i<s;++i)
		a[i] -= b[i];

}

// print vector using cout
void print_vec(std::vector<double>& v,int N = 0)
{
	if (N==0) N = v.size();
	for(int i=0;i<N;++i)
		std::cout << v[i] << std::endl;
}

// write vector to a file named fileName.
// writes values separated by newline
// write N values, default N = v.size()
void write_vec(std::vector<double> v,std::string fileName, int N = 0)
{
	if(N==0) N = v.size();
	std::ofstream out;
	out.open(fileName);
	for(int i=0;i<N;++i)
		out<<v[i] << '\n';
}


// write vector to a file named fileName.
// writes values separated by newline
// write N values, default N = v.size()
void write_vec(std::vector<int> v,std::string fileName, int N = 0)
{
	if(N==0) N = v.size();
	std::ofstream out;
	out.open(fileName);
	for(int i=0;i<N;++i)
		out<<v[i] << '\n';
}



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



void rand_vecs(std::vector<std::vector<double> >& r,
	int N, int dim, double min, double max,
	std::default_random_engine& gen, double norm = -1.)
{
	std::uniform_real_distribution<double> dist(min,max);
	for(int i=0;i<N;++i) {
		for(int d=0;d<dim;++d) {
			r[i][d] = dist(gen);
		}
		if(norm > 0.)
			normalize(r[i],norm);
	}
}



#endif
