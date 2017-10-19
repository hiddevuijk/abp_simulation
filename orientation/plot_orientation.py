import numpy as np
import matplotlib.pyplot as plt
from sys import exit

bins = np.loadtxt("pAvg_bins.dat")



plt.subplot(3,3,1)
plt.title(r"$p_x(x)$")
pz = np.loadtxt("pAvgx.dat",delimiter=';')
pzz = pz[:,0]
plt.plot(bins,pzz)
plt.grid(True)

plt.subplot(3,3,2)
plt.title(r"$p_y(x)$")
pz = np.loadtxt("pAvgx.dat",delimiter=';')
pzz = pz[:,1]
plt.plot(bins,pzz)
plt.grid(True)

plt.subplot(3,3,3)
plt.title(r"$p_z(x)$")
pz = np.loadtxt("pAvgx.dat",delimiter=';')
pzz = pz[:,2]
plt.plot(bins,pzz)
plt.grid(True)

plt.subplot(3,3,4)
plt.title(r"$p_x(y)$")
pz = np.loadtxt("pAvgy.dat",delimiter=';')
pzz = pz[:,0]
plt.plot(bins,pzz)
plt.grid(True)

plt.subplot(3,3,5)
plt.title(r"$p_y(y)$")
pz = np.loadtxt("pAvgy.dat",delimiter=';')
pzz = pz[:,1]
plt.plot(bins,pzz)
plt.grid(True)

plt.subplot(3,3,6)
plt.title(r"$p_z(y)$")
pz = np.loadtxt("pAvgy.dat",delimiter=';')
pzz = pz[:,2]
plt.plot(bins,pzz)
plt.grid(True)


plt.subplot(3,3,7)
plt.title(r"$p_x(z)$")
pz = np.loadtxt("pAvgz.dat",delimiter=';')
pzz = pz[:,0]
plt.plot(bins,pzz)
plt.grid(True)


plt.subplot(3,3,8)
plt.title(r"$p_y(z)$")
pz = np.loadtxt("pAvgz.dat",delimiter=';')
pzz = pz[:,1]
plt.plot(bins,pzz)
plt.grid(True)

plt.subplot(3,3,9)
plt.title(r"$p_z(z)$")
pz = np.loadtxt("pAvgz.dat",delimiter=';')
pzz = pz[:,2]
plt.plot(bins,pzz)
plt.grid(True)



plt.tight_layout()
plt.show()






