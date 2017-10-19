import numpy as np

v0 = 0.
navg = 5

rm_vh = 0.5
dr_vh = 0.2

Nr_vh = rm_vh/dr_vh + 1


N = 102;
Dt = 1./1;
Dr = 20.*Dt;
gamma = 1.;
beta = 1.;
eps = 1.;
sigma = 1.;
seed = 123456789;
d = sigma
dt = (1.e-4)*d*d/Dt;

tm_vh = .5*d*d/Dt
dt_vh = .05*d*d/Dt;
Nt_vh = tm_vh/dt_vh + 1;

teq = 0.*d*d/Dt;

rho = .1

L = (N/rho)**(1/3.)


name = ""
infile = open("input.txt",'w')
infile.write("N=\n%i\n" % N)
infile.write("rho=\n%f\n" % rho)
infile.write("Dt=\n%f\n" % Dt)
infile.write("Dr=\n%f\n" % Dr)
infile.write("v0=\n%f\n" % v0)
infile.write("gamma=\n%f\n" % gamma)
infile.write("beta=\n%f\n" % beta)
infile.write("eps=\n%f\n" % eps)
infile.write("sigma=\n%f\n" % sigma)
infile.write("dt=\n%f\n" % dt)
infile.write("dt_vh=\n%f\n" % dt_vh)
infile.write("tm_vh=\n%f\n" % tm_vh)
infile.write("dr_vh=\n%f\n" % dr_vh)
infile.write("rm_vh=\n%f\n" % rm_vh)
infile.write("teq=\n%i\n" % teq)
infile.write("navg=\n%i\n" % navg)
infile.write("seed=\n%i\n" % seed)
infile.write("name=\n" + name)

infile.close()






