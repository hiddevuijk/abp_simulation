import numpy as np
import matplotlib.pyplot as plt
from sys import exit


bins = np.loadtxt("bins.dat")[:-1]
bs = bins[1] - bins[0]

rhox = np.loadtxt("rhox.dat")[:-1]
mx = np.mean(rhox)
rhox = rhox/mx - np.ones(rhox.shape)

rhoy = np.loadtxt("rhoy.dat")[:-1]
my = np.mean(rhoy)
rhoy = rhoy/my - np.ones(rhoy.shape)
rhoz = np.loadtxt("rhoz.dat")[:-1]
mz = np.mean(rhoz)
rhoz = rhoz/mz - np.ones(rhoz.shape)



plt.plot(bins,rhoy)

plt.xlabel('y')
plt.tight_layout()
plt.show()






