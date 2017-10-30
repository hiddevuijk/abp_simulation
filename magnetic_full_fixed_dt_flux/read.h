#ifndef GUARD_read_h
#define GUARD_read_h

#include <string>
#include <fstream>
#include <iostream>

using namespace std;

void read_variables_vhove(int& N,double& rho, double& Dt,double& Dr,double& v0, double& gamma,
	double& beta,double& eps,double& sigma,double& dt,double& dt_vh, double& tm,
	double& dr_vh, double& rm,double& teq,int& navg, int& seed,
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
	infile >> v0;
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
	infile >> dt_vh;
	infile >>temp;
	infile >> tm;
	infile >> temp;
	infile >> dr_vh;
	infile >> temp;
	infile >> rm;
	infile >> temp;
	infile >> teq;
	infile >> temp;
	infile >> navg;
	infile >> temp;
	infile >> seed;
	infile >>temp;
	infile >> name;
	
}

void read_variables_pairdist(int& N,double& rho, double& Dt,double& Dr,double& v0,
	double& gamma,double& beta,double& eps,double& sigma,double& dt,
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
	infile >> Dr;
	infile >>temp;
	infile >> v0;
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
	infile >> bs;
	infile >> temp;
	infile >> seed;
	infile >>temp;
	infile >> name;
	
}


void read_variables_temp(int& N,double& rho, double& Dt,double& Dr,double& v0,double& Tbase,
	double& T,double& omega0,double& beta,double& eps,double& sigma,double& dt,
	double& tf,double& teq, int& navg,double& bs_pAvg, int& seed,
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
	infile >> v0;
	infile >>temp;
	infile >> Tbase;
	infile >>temp;
	infile >> T;
	infile >>temp;
	infile >> omega0;
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
	infile >> bs_pAvg;
	infile >> temp;
	infile >> seed;
	infile >>temp;
	infile >> name;
	
}

void read_variables_magnetic(int& N,double& rho, double& Dt,double& Dr,double& v0,
	double& m, double& qB, double& omega0,double& beta,double& eps,
	double& sigma,double& dt, double& tf,double& teq, int& navg,double& bs_pAvg,
	int& seed,std::string& name,std::string fileName)
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
	infile >> v0;
	infile >>temp;
	infile >> m;
	infile >>temp;
	infile >> qB;
	infile >>temp;
	infile >> omega0;
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
	infile >> bs_pAvg;
	infile >> temp;
	infile >> seed;
	infile >>temp;
	infile >> name;
	
}

void read_variables_magnetic(int& N,double& rho, double& Dt,double& Dr,double& v0,double& qB,
	double& omega0,double& beta,double& eps,double& sigma,double& dt,
	double& tf,double& teq, int& navg,double& bs_pAvg, int& seed,
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
	infile >> v0;
	infile >>temp;
	infile >> qB;
	infile >>temp;
	infile >> omega0;
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
	infile >> bs_pAvg;
	infile >> temp;
	infile >> seed;
	infile >>temp;
	infile >> name;
	
}



void read_variables_magnetic_vanHove(int& N,double& rho, double& Dt,double& Dr,double& v0,double& qB,
	double& omega0,double& beta,double& eps,double& sigma,double& dt,
	double& tf,double& teq, int& navg,double& bs,double& tbs, double& tmax, int& seed,
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
	infile >> v0;
	infile >>temp;
	infile >> qB;
	infile >>temp;
	infile >> omega0;
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
	infile >> bs;
	infile >> temp;
	infile >> tbs;
	infile >> temp;
	infile >> tmax;
	infile >> temp;
	infile >> seed;
	infile >>temp;
	infile >> name;
	
}

void read_variables_orientation(int& N,double& rho, double& Dt,double& Dr,double& v0,
	double& omega0, double& gamma,double& beta,double& eps,double& sigma,double& dt,
	double& tf,double& teq, int& navg,double& bs_pAvg, int& seed,
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
	infile >> v0;
	infile >>temp;
	infile >> omega0;
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
	infile >> bs_pAvg;
	infile >> temp;
	infile >> seed;
	infile >>temp;
	infile >> name;
	
}


void read_variables_diffconst(int& N,double& rho, double& Dt,double& Dr,double& v0,
	double& gamma,double& beta,double& eps,double& sigma,double& dt,
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
	infile >> v0;
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

