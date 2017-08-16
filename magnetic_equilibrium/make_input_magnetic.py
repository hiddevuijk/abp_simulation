import numpy as np

omega0 = 2.
qB = 1.
navg = 1000;
bs =0.05;
N = 500;
Dt = 1./1;
gamma = 1.;
seed = 123456789;
d = 1.;
dt = (5.e-4)*d*d/Dt;
tf = .5*d*d/Dt;
teq = 5.*d*d/Dt;

rho = .5

L = (N/rho)**(1/3.)


name = ""
infile = open("input.txt",'w')
infile.write("N=\n%i\n" % N)
infile.write("rho=\n%f\n" % rho)
infile.write("Dt=\n%f\n" % Dt)
infile.write("qB=\n%f\n" % qB)
infile.write("omega0\n%f\n" % omega0)
infile.write("gamma=\n%f\n" % gamma)
infile.write("dt=\n%f\n" % dt)
infile.write("tf=\n%f\n" % tf)
infile.write("teq=\n%i\n" % teq)
infile.write("navg=\n%i\n" % navg)
infile.write("bs=\n%f\n" % bs)
infile.write("seed=\n%i\n" % seed)
infile.write("name=\n" + name)

infile.close()






