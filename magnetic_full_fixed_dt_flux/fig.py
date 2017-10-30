import numpy as np
import matplotlib.pyplot as plt
from sys import exit

b = np.loadtxt("bins.dat")
bs = b[1] - b[0]
rx0 = np.loadtxt("rhox.dat")/bs
ry0 = np.loadtxt("rhoy.dat")/bs
rz0 = np.loadtxt("rhoz.dat")/bs

#rx0 = np.loadtxt("fluxy.dat")/bs
#ry0 = np.loadtxt("fluyx.dat")/bs
#rz0 = np.loadtxt("fluyy.dat")/bs



plt.subplot(3,1,1)
plt.plot(b[:-1],rx0[:-1])


plt.subplot(3,1,2)
plt.plot(b[:-1],ry0[:-1])
plt.subplot(3,1,3)
plt.plot(b[:-1],rz0[:-1])


plt.tight_layout()
plt.show()






