import numpy as np

v0 = 5.
omega0 = 1.
qB = 2.
navg = 10;
nbin = 100;
N = 30;
Dt = 1./1;
Dr = 20.*Dt;
beta = 1.;
eps = 0.;
sigma = 1.;
seed = 123456789;
d = sigma
dt = (1.e-3)*d*d/Dt;
tf = .5*d*d/Dt;
teq = 10.*d*d/Dt;

rho = 0.2

L = (N/rho)**(1/3.)

bs_pAvg =L/nbin;

name = ""
infile = open("input.txt",'w')
infile.write("N=\n%i\n" % N)
infile.write("rho=\n%f\n" % rho)
infile.write("Dt=\n%f\n" % Dt)
infile.write("Dr=\n%f\n" % Dr)
infile.write("v0=\n%f\n" % v0)
infile.write("qB=\n%f\n" % qB)
infile.write("omega0\n%f\n" % omega0)
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






