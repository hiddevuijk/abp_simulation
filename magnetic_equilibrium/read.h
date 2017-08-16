#ifndef GUARD_read_h
#define GUARD_read_h

#include <string>
#include <fstream>
#include <iostream>

using namespace std;



void read_variables_magnetic(int& N,double& rho, double& Dt,double& qB,
	double& omega0, double& gamma,double& dt,
	double& tf,double& teq, int& navg,double& bs, int& seed,
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
	infile >> qB;
	infile >>temp;
	infile >> omega0;
	infile >>temp;
	infile >> gamma;
	infile >>temp;
	infile >> dt;
	infile >>temp;
	infile >> tf;
	infile >>temp;
	infile >> teq;
	infile >> temp;
	infile >> navg;
	infile >> temp;
	infile >> bs;
	infile >> temp;
	infile >> seed;
	infile >>temp;
	infile >> name;
	
}

#endif

