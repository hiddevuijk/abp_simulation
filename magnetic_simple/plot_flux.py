import numpy as np
import matplotlib.pyplot as plt
from sys import exit

Dr = 20


fxyX = np.loadtxt("fxyX.dat",delimiter=';')
fxyY = np.loadtxt("fxyY.dat",delimiter=';')

fx = np.mean(fxyX,axis=0)
fy = np.mean(fxyY,axis=0)
N = fxyX.shape[0]

x,y = np.meshgrid(np.linspace(0,10,N),np.linspace(0,10,N))

plt.subplot(2,1,1)
plt.quiver(x,y,fxyX,fxyY,units="width")
#qk = plt.quiverkey(Q,0.9,0.9,2.,"label",labelpos='E',coordinates='figure')

plt.subplot(2,1,2)
plt.plot(np.linspace(0,10.,fx.shape[0]),fx,label='flux X')
plt.plot(np.linspace(0,10.,fy.shape[0]),fy,label='flux Y')
plt.legend()

plt.tight_layout()
plt.show()






