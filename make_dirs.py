import numpy as np


N = 100;
Dt = 1./1;
Dr = 1./10.;
gamma = 1.;
beta = 1.;
eps = 1.;
sigma = 1.;
L = 10.;
seed = 123456789;
dt = (1.e-5)*sigma*sigma/Dt;
tf = 1e1;
teq = 1.;


name = ""
infile = open("input.txt",'w')
infile.write("N=\n%i\n" % N)
infile.write("Dt=\n%f\n" % Dt)
infile.write("Dr=\n%f\n" % Dr)
infile.write("gamma=\n%f\n" % gamma)
infile.write("beta=\n%f\n" % beta)
infile.write("eps=\n%f\n" % eps)
infile.write("sigma=\n%f\n" % sigma)
infile.write("L=\n%f\n" % L)
infile.write("dt=\n%f\n" % dt)
infile.write("tf=\n%f\n" % tf)
infile.write("teq=\n%i\n" % teq)
infile.write("sees=\n%i\n" % seed)
infile.write("name=\n" + name)

infile.close()






