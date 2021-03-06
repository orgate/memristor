# From calculation, we expect that the local minimum occurs at x=9/4

from numpy import *
import matplotlib.pyplot as plt

#x_old = 0
#x_new = 2 # The algorithm starts at x=6
eps = 0.001 # step size
precision = 0.0000001

xarray = arange(0,1,eps)
phin = xarray*0
phin[-1] = 1.0
#def f_prime(x):
#    return 4 * x**3 - 9 * x**2

#in forward direction
def f(prev,next):
	value = 0.026*log(abs(80*(next-prev)/(-2*100*eps)))
#	value = 100.0*eps/(80.0*exp(phi/0.026)) + prev
#	value = -exp(x)
	return value
 


for i in range(len(xarray)):
	if (i>0)&(i<len(xarray)-1):
		prev = phin[-i-2]
		next = phin[-i]
		phin[-i-1] = f(prev,next)
		print "phin[-i-1] is ",phin[-i-1]

for i in range(len(xarray)):
	if (i>0)&(i<len(xarray)-1):
		prev = phin[i-1]
		next = phin[i+1]
#		phi_old = phin[i]
#		phi_new = 1.0	

#		tempphi = 
#		phi_old = phi_new
#	x_new = f(x_old)
#		phi_new = tempphi
		phin[i] = f(prev,next)
		print "phin[i] is ",phin[i]


#		print "phi_new is ",phi_new

#print "Local minimum occurs at ", phi_new
#plt.plot(xarray,phin)
plt.show()
