import numpy as np
import matplotlib.pyplot as plt
from sys import exit


fxyX = np.loadtxt("fxyX.dat",delimiter=';')
fxyY = np.loadtxt("fxyY.dat",delimiter=';')

N = fxyX.shape[0]

x,y = np.meshgrid(np.linspace(0,10,N),np.linspace(0,10,N))

plt.figure()
Q = plt.quiver(x,y,fxyX,fxyY,units="width")
qk = plt.quiverkey(Q,0.9,0.9,2.,"label",labelpos='E',coordinates='figure')
plt.show()






