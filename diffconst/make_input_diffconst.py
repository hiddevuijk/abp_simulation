import numpy as np

rho = .7
v0 = 0.
N = 100;
T = 100; # end time in tau_b
Nt = 100;
Dt = 1./1;
Dr = 20.*Dt;
gamma = 1.;
beta = 1.;
eps = 1.;
sigma = 1.;
seed = 123456780;
d = sigma
dt = (1.e-5)*d*d/Dt;
teq = 10.*d*d/Dt;



tf = 1.*T/Nt
tf *= d*d/Dt;

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
infile.write("tf=\n%f\n" % tf)
infile.write("teq=\n%i\n" % teq)
infile.write("Nt=\n%i\n" % Nt)
infile.write("seed=\n%i\n" % seed)
infile.write("name=\n" + name)

infile.close()






