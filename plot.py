import numpy as np
import matplotlib.pyplot as plt
from sys import exit



rhox = np.loadtxt("rhox.dat")[:-1]
rhoy = np.loadtxt("rhoy.dat")[:-1]
rhoz = np.loadtxt("rhoz.dat")[:-1]

bins = np.loadtxt("bins.dat")[:-1]

plt.subplot(3,1,1)
plt.plot(bins,rhox)

plt.subplot(3,1,2)
plt.plot(bins,rhoy)

plt.subplot(3,1,3)
plt.plot(bins,rhoz)


plt.tight_layout()
plt.show()






