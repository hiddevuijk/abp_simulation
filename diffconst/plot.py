import numpy as np
import matplotlib.pyplot as plt
from sys import exit



inp = open("input.txt",'r')
params = inp.readlines()
N = int(params[1])
rho = float(params[3])
Dt = float(params[5])
Dr = float(params[7])
gamma = float(params[9])
beta = float(params[11])
eps = float(params[13])
sigma = float(params[15])
dt = float(params[17])
tf = float(params[19])
teq = float(params[21])
Nt = float(params[23])
seed = float(params[25])
name = params[27]

T = Nt*tf

dr = np.loadtxt(name+".dat")
t = np.arange(tf,T+tf,tf)
theory = np.arange(tf,T+tf,tf)*6*Dt

plt.plot(t,dr,label="sim")
plt.plot(t,theory,label="6*Dt*t")
plt.legend()
plt.tight_layout()
plt.show()

