import numpy as np

v0 = 0.
omega0 = 1.
qB = 10.
navg = 1000;
bs_pAvg =0.1;
N = 100;
Dt = 1./1;
Dr = 20.*Dt;
gamma = 1.;
beta = 1.;
eps = 0.;
sigma = 1.;
seed = 123456701;
d = sigma
dt = (1.e-3)*d*d/Dt;
tf = .5*d*d/Dt;
teq = 10.*d*d/Dt;

rho = .1

L = (N/rho)**(1/3.)


name = ""
infile = open("input.txt",'w')
infile.write("N=\n%i\n" % N)
infile.write("rho=\n%f\n" % rho)
infile.write("Dt=\n%f\n" % Dt)
infile.write("Dr=\n%f\n" % Dr)
infile.write("v0=\n%f\n" % v0)
infile.write("qB=\n%f\n" % qB)
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






