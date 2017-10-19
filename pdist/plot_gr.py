import numpy as np
import matplotlib.pyplot as plt
from sys import exit



g = np.loadtxt("g.dat")
r = np.loadtxt("gr.dat")

plt.plot(r,g)
plt.show()





