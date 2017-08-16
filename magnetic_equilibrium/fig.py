import numpy as np
import matplotlib.pyplot as plt
from sys import exit

b = np.loadtxt("bins.dat")
bs = b[1] - b[0]
rx0 = np.loadtxt("rhox.dat")/bs
ry0 = np.loadtxt("rhoy.dat")/bs
rz0 = np.loadtxt("rhoz.dat")/bs

#rx = np.loadtxt("rhox_.dat")/bs
#ry = np.loadtxt("rhoy_.dat")/bs
#rz = np.loadtxt("rhoz_.dat")/bs

plt.plot(b,ry0,label="f1")
#plt.plot(b,ry,label="S as I")
plt.title("y")
#plt.ylim([0.035,0.06])
plt.legend()
plt.tight_layout()
plt.show()






