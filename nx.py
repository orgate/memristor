# From calculation, we expect that the local minimum occurs at x=9/4

from numpy import *
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.special import gamma, airy
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

#Y0 = -exp(-5.0334620275)
Y0 = -exp(-5.033462027062)
Y1 = 1.0
Y = [Y0, Y1]
k=-0.026

def func(y, t):
	return [y[0]*y[0]/k,y[0]]

#def func(y, t):
#	return [1*exp(y[0]/k)/80.0,y[0]]


eps = 0.01 # step size
x = arange(0, 4.0, eps)
xarray = arange(0, 4.0, eps)
size = len(x)
phin = odeint(func, Y, x)
anscheck = xarray*0
phin1 = xarray*0
#phin[-1] = 0.04

print phin[-1][1]

for i in range(len(xarray)):
	phin1[i]=phin[i][1]

#for i in range(len(xarray)):
#	if (i!=0):
#		anscheck[i] = phin[i][0]*exp(phin[i][1]/0.026)
#nconc = xarray*0
#for i in range(len(xarray)):
#	nconc[i] = exp(phin[i][1]/0.026)
#print anscheck

#plt.title('Values of electron concn.(green) and phin(blue)')
#plt.title('Phin')
#plt.title('n(x)')
#plt.xlabel('length x/L')
#plt.ylabel('n (in per cubic cm) and electron fermi potential (in V)')
#plt.ylabel('n(x)')
#plt.ylabel('phin')
#plt.plot(xarray,phin1)
#plt.plot(xarray,anscheck)
#plt.plot(xarray,nconc)
#plt.show()
