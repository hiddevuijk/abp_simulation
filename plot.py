import numpy as np
import matplotlib.pyplot as plt
from sys import exit

binsize = 0.05


infile = open("input.txt")
input_vals = infile.readlines()
N = int(input_vals[1][:-1])
rho = float(input_vals[3][:-1])
L = (N/rho)**(1./3)


pd_mat = np.loadtxt("pdist_mat.dat",delimiter=';')
nv = pd_mat.shape[0]

pdist_r = np.arange(binsize,L,binsize)
g = np.zeros(pdist_r.size)

for nvi in range(nv):
	gtemp = np.zeros(pdist_r.size)
	pdistances = pd_mat[nvi]	
	di = 0

	for i in range(len(pdist_r)):
		if di+1 < pdistances.shape[0]:
			while pdistances[di] < pdist_r[i]:
				gtemp[i] += 1.
				di += 1
				if di>=pdistances.shape[0]:
					break
		else:
			gtemp[i] += 1.

	for i in range(len(pdist_r)):
		if nvi ==0 :pdist_r[i] -= 0.5*binsize
		gtemp[i] /= rho
		print pdist_r[i]
	
		gtemp[i] /=4*np.pi*pdist_r[i]*pdist_r[i]*binsize
		gtemp[i] /= 0.5*(N-1)
	g += gtemp/pd_mat.shape[0]
plt.plot(pdist_r,g)
plt.show()




'''
intg = 0
for i in range(len(g)):
	intg += g[i]*rho*4*np.pi*pdist_r[i]*pdist_r[i]*binsize
print intg
print N*(N-1)/2
'''



'''
pdistances = np.loadtxt('pdist.dat')

pdist_r = np.arange(binsize,L,binsize)
g = np.zeros(pdist_r.size)
di = 0

for i in range(len(pdist_r)):
	if di+1 < pdistances.size:
		while pdistances[di] < pdist_r[i]:
			g[i] += 1.
			di += 1
	else:
		g[i] += 1.

for i in range(len(pdist_r)):
	pdist_r[i] -= 0.5*binsize
	g[i] /= rho
	g[i] /=4*np.pi*pdist_r[i]*pdist_r[i]*binsize
	g[i] /= 0.5*(N-1)
'''
intg = 0
for i in range(len(g)):
	intg += g[i]*rho*4*np.pi*pdist_r[i]*pdist_r[i]*binsize
print intg
print N*(N-1)/2
'''


plt.axvline(L)
plt.plot(pdist_r,g)
plt.show()

'''


