#ifndef GUARD_read_h
#define GUARD_read_h

#include <string>
#include <fstream>
#include <iostream>

using namespace std;


void read_variables(int& N,double& Dt,double& Dr,double& gamma,
	double& beta,double& eps,double& sigma,double& L,double& dt,
	double& tf,double& teq,int& seed,std::string& name,
	std::string fileName)
{
	std::string temp;
	std::ifstream infile(fileName);

	infile >>temp;
	infile >> N;
	infile >>temp;
	infile >> Dt;
	infile >>temp;
	infile >> Dr;
	infile >>temp;
	infile >> gamma;
	infile >>temp;
	infile >> beta;
	infile >>temp;
	infile >> eps;
	infile >>temp;
	infile >> sigma;
	infile >>temp;
	infile >> L;
	infile >>temp;
	infile >> dt;
	infile >>temp;
	infile >> tf;
	infile >>temp;
	infile >> teq;
	infile >>temp;
	infile >> seed;
	infile >>temp;
	infile >> name;
	
}


#endif

