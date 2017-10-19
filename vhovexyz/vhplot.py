import numpy as np
import matplotlib.pyplot as plt
from sys import exit

vh = np.loadtxt("vhove.dat",delimiter=';')[0]
r = np.loadtxt("r.dat")
print 4*np.pi*(r[1]-r[0])*np.dot(vh,np.multiply(r,r))
plt.plot(r,vh)
plt.show()


