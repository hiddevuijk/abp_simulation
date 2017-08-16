import numpy as np
import matplotlib.pyplot as plt
from sys import exit

q = 1.
rho = .5
N = 5000

omega0 = 1.

L = (N/rho)**(1./3)
omega = omega0*2*np.pi/L
n = 5000
gamma = 1.

def w(y):
	return q*np.sin(omega*y)

def wp(y):
	return omega*q*np.cos(omega*y)


def a(y):
	return w(y)/(gamma*gamma+w(y)*w(y))

def b(y):
	return gamma/(gamma*gamma+w(y)*w(y))

def A(y):
	x = gamma*gamma*wp(y)-wp(y)*w(y)*w(y)
	return x/((gamma*gamma+w(y)*w(y))**2)

def B(y):
	x = -2.*gamma*wp(y)*w(y)
	return x/((gamma*gamma+w(y)*w(y))**2)

def Ay(y):
	return (b(y)*B(y)+a(y)*A(y))


Y = np.linspace(0,L,n)
plt.plot(Y,Ay(Y),label="+")
omega0=-1
plt.plot(Y,Ay(Y),label="-")

plt.plot(Y,np.zeros(len(Y)))

plt.legend()
plt.tight_layout()
plt.show()



