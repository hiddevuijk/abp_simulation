import numpy as np

navg = 20;
N = 100;
Dt = 1./1;
Dr = 20.*Dt;
gamma = 1.;
beta = 1.;
eps = 1.;
sigma = 1.;
seed = 123456780;
d = sigma
dt = (1.e-5)*d*d/Dt;
tf = 10.*d*d/Dt;
teq = 10.*d*d/Dt;

rho = 0.2

L = (N/rho)**(1/3.)

name = ""
infile = open("input.txt",'w')
infile.write("N=\n%i\n" % N)
infile.write("rho=\n%f\n" % rho)
infile.write("Dt=\n%f\n" % Dt)
infile.write("Dr=\n%f\n" % Dr)
infile.write("gamma=\n%f\n" % gamma)
infile.write("beta=\n%f\n" % beta)
infile.write("eps=\n%f\n" % eps)
infile.write("sigma=\n%f\n" % sigma)
infile.write("dt=\n%f\n" % dt)
infile.write("tf=\n%f\n" % tf)
infile.write("teq=\n%i\n" % teq)
infile.write("navg=\n%i\n" % navg)
infile.write("seed=\n%i\n" % seed)
infile.write("name=\n" + name)

infile.close()






