# From calculation, we expect that the local minimum occurs at x=9/4

from numpy import *
import matplotlib.pyplot as plt

#x_old = 0
#x_new = 2 # The algorithm starts at x=6
eps = 0.001 # step size
precision = 0.0000001

xarray = arange(0,1,eps)
phin = xarray*0
#phin[-1] = 1.0
#def f_prime(x):
#    return 4 * x**3 - 9 * x**2

#in forward direction
def f(phi,prev):
	value = 100.0*eps/(80.0*exp(phi/0.026)) + prev
#	value = -exp(x)
	return value
 
#in backward direction
def f1(phi):
	value = -100.0*eps/(80.0*exp(phi/0.026)) + phi
#	value = -exp(x)
	return value


for i in range(len(xarray)):
	if (i>0)&(i<len(xarray)-1):
		prev = phin[i-1]
		phi_old = phin[i]
		phi_new = 1.0	

		while abs(phi_new - phi_old) > precision:
			tempphi = f(phi_old,prev)
			phi_old = phi_new
#	x_new = f(x_old)
			phi_new = tempphi
		phin[i] = phi_new
plt.plot(xarray,phin)

for i in range(len(xarray)):
	if (i>0)&(i<len(xarray)-1):
		prev = phin[-i-1]
		phi = phin[-i]
		tempphiprev = f1(phi)
		phin[-i-1] = tempphiprev


#		print "phi_new is ",phi_new

print "Local minimum occurs at ", phi_new
#plt.plot(xarray,phin)
plt.show()
