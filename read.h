#ifndef GUARD_read_h
#define GUARD_read_h

#include <string>
#include <fstream>
#include <iostream>

using namespace std;


void read_variables_pairdist(int& N,double& rho, double& Dt,double& Dr,double& gamma,
	double& beta,double& eps,double& sigma,double& dt,
	double& tf,double& teq, int& navg, int& seed,
	std::string& name,std::string fileName)
{
	std::string temp;
	std::ifstream infile(fileName);

	infile >>temp;
	infile >> N;
	infile >> temp;
	infile >> rho;
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
	infile >> dt;
	infile >>temp;
	infile >> tf;
	infile >>temp;
	infile >> teq;
	infile >> temp;
	infile >> navg;
	infile >> temp;
	infile >> seed;
	infile >>temp;
	infile >> name;
	
}



void read_variables_diffconst(int& N,double& rho, double& Dt,double& Dr,double& gamma,
	double& beta,double& eps,double& sigma,double& dt,
	double& tf,double& teq, int& Nt, int& seed,
	std::string& name,std::string fileName)
{
	std::string temp;
	std::ifstream infile(fileName);

	infile >>temp;
	infile >> N;
	infile >> temp;
	infile >> rho;
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
	infile >> dt;
	infile >>temp;
	infile >> tf;
	infile >>temp;
	infile >> teq;
	infile >> temp;
	infile >> Nt;
	infile >> temp;
	infile >> seed;
	infile >>temp;
	infile >> name;
	
}


#endif

