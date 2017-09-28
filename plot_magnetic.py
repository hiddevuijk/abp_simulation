import numpy as np
import matplotlib.pyplot as plt
from sys import exit


q = 1
omega0 = 1
rho = 0.2
N = 300

L = (N/rho)**(1./3)
omega = 2*np.pi*omega0/L

def w(y):
	return q*np.sin(omega*y)
def wp(y):
	return q*omega*np.cos(omega*y)

def Ap(y):
	D = 1+w(y)*w(y)
	return wp(y)*(1-w(y)*w(y))/(D*D)
def Bp(y):
	D= 1+w(y)*w(y)
	return -2.*w(y)*wp(y)/(D*D)

y = np.linspace(0,L,75)
px = -1.*Ap(y)
py = -1.*Bp(y)

plt.plot(y,px/12,label='px')
plt.plot(y,py/12,label='py')

plt.legend()
plt.tight_layout()

plt.show()




