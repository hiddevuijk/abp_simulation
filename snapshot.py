import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sys import exit


def drawSphere(xc,yc,zc,r):
	u,v = np.mgrid[0:2*np.pi:30j,0:np.pi:30j]
	x = np.cos(u)*np.sin(v)
	y = np.sin(u)*np.sin(v)
	z = np.cos(v)
	x = r*x+xc
	y = y*r+yc
	z = z*r+zc
	return (x,y,z)

def drawBox(L,ax):
	ax.plot_wireframe([0,0],[0,L],[0,0])

def snapshot(filename,sigma,L):
	data = np.loadtxt(filename,delimiter=';')
	x,y,z = data[:,0],data[:,1],data[:,2]
	fig = plt.figure()
	ax = fig.add_subplot(111,projection='3d')

	for (xi,yi,zi) in zip(x,y,z):
		(xs,ys,zs) = drawSphere(xi,yi,zi,sigma)
		ax.plot_wireframe(xs,ys,zs,color=np.random.rand(3,))
		#ax.plot_surface(xs,ys,zs,color=np.random.rand(3,))

	#drawBox(L,ax)

	ax.set_xlabel("x")
	ax.set_ylabel("y")
	ax.set_zlabel("z")

	ax.set_xlim([0,L])
	ax.set_ylim([0,L])
	ax.set_zlim([0,L])
	plt.title("N=100")
	plt.show()

N = 100
rho = 0.7
L = (N/rho)**(1/3.)

snapshot('r100.dat',1.,L)




