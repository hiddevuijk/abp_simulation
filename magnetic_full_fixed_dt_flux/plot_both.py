import numpy as np
import matplotlib.pyplot as plt
from sys import exit

v0 = 10.
Dr = 20.
q = 5.

bins = np.loadtxt("bins.dat")[:-1]
bs = (bins[1] - bins[0])



pyx = np.loadtxt("pAvgy.dat",delimiter=';')[:,0][:-1]*bs
pyy = np.loadtxt("pAvgy.dat",delimiter=';')[:,1][:-1]*bs
plt.scatter(bins,pyx,label=r'$p_x(y)$',color='blue')
plt.scatter(bins,pyy,label=r'$p_y(y)$',color='red')


L = 10
omega = 2*np.pi/L


def w(y,q,omega):
	return q*np.sin(omega*y)
def wp(y,q,omega):
	return q*omega*np.cos(omega*y)

def Ap(y,q,omega):
	D = 1+w(y,q,omega)**2
	return wp(y,q,omega)*(1-w(y,q,omega)**2)/(D*D)
def Bp(y,q,omega):
	D= 1+w(y,q,omega)**2
	return -2.*w(y,q,omega)*wp(y,q,omega)/(D*D)
def px(y,v0,Dr,q,omega):
	return -v0*Ap(y,q,omega)/(6*Dr)
def py(y,v0,Dr,q,omega):
	return -v0*Bp(y,q,omega)/(6*Dr)

y = np.linspace(0,L,1000)
px_y = px(y,v0,Dr,q,omega)
py_y = py(y,v0,Dr,q,omega)

plt.plot(y,px_y,color='blue')
plt.plot(y,py_y,color='red')


ma = max(max(pyx),max(pyy),max(px_y),max(py_y))
mi = min(min(pyx),min(pyy),min(px_y),min(py_y))
ma = max(ma,-1*mi)
m = 1.1*ma

plt.ylim([-m,m])
plt.xlim([0,L])
plt.xlabel('y')
plt.legend()
plt.tight_layout()
plt.show()






