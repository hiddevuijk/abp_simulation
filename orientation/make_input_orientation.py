import numpy as np
from sys import exit

v0 = 20.;
omega0 = 1.
navg = 50;
bs_pAvg =.5;



N = 500;
Dt = 1./1;
Dr = 25.*Dt;
gamma = 1.;
beta = 1.;
eps = 1.;
sigma = 1.;
seed = 123456789;

d = sigma
dt = (1.e-3)*d*d/Dt;
tf = 1.*d*d/Dt;
teq = 5.*d*d/Dt;

rho = .2
L = (N/rho)**(1/3.)
#L = 50.
#rho = N/(L**3.)


name = ""
infile = open("input.txt",'w')
infile.write("N=\n%i\n" % N)
infile.write("rho=\n%f\n" % rho)
infile.write("Dt=\n%f\n" % Dt)
infile.write("Dr=\n%f\n" % Dr)
infile.write("v0=\n%f\n" % v0)
infile.write("omega0\n%f\n" % omega0)
infile.write("gamma=\n%f\n" % gamma)
infile.write("beta=\n%f\n" % beta)
infile.write("eps=\n%f\n" % eps)
infile.write("sigma=\n%f\n" % sigma)
infile.write("dt=\n%f\n" % dt)
infile.write("tf=\n%f\n" % tf)
infile.write("teq=\n%i\n" % teq)
infile.write("navg=\n%i\n" % navg)
infile.write("bs_pAvg=\n%f\n" % bs_pAvg)
infile.write("seed=\n%i\n" % seed)
infile.write("name=\n" + name)

infile.close()






