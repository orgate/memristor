# From calculation, we expect that the local minimum occurs at x=9/4

from numpy import *
import matplotlib.pyplot as plt

#x_old = 0
#x_new = 2 # The algorithm starts at x=6
eps = 0.01 # step size
precision = 0.00001

xarray = arange(0,1,eps)
phin = xarray*0
anscheck = xarray*0
phin[-1] = 0.04
#def f_prime(x):
#    return 4 * x**3 - 9 * x**2

##def f(phi,prev,next):
##	value = ( prev+next*( exp(( next-phi )/0.026 )))/( 1+exp(( next-phi )/0.026 ) )
#	value = -exp(x)
##	return value
 
def f(phi,prev,next,precision):
#	value = ( prev+next*( exp(( next-phi )/0.026 )))/( 1+exp(( next-phi )/0.026 ) )
	value = ( prev+next*( exp(( next-prev )/0.026 )))/( 1+exp(( next-prev )/0.026 ) )
	if abs(phi-value)>precision:
		f(value,prev,next,precision)
	return value

for m in range(100):
	for i in range(len(xarray)):
		if (i>0)&(i<len(xarray)-1):
			prev = phin[i-1]
			next = phin[i+1]
			phi_old = phin[i]
##			phi_new = 10	

			iterate = 0
##			while (abs(phi_new - phi_old) > precision)&(iterate<10):
##				tempphi = f(phi_old,prev,next)
##				phi_old = phi_new
#	x_new = f(x_old)
##				phi_new = tempphi
##				iterate+=1
##			phin[i] = phi_new
			phin[i] = f(phi_old,prev,next,precision)

#print phin

	for i in range(len(xarray)):
		if (i>0)&(i<len(xarray)-1):
			prev = phin[-i-2]
			next = phin[-i]
			phi_old = phin[-i-1]
##			phi_new = 10	

##			while (abs(phi_new - phi_old) > precision)&(iterate<10):
##				tempphi = f(phi_old,prev,next)
##				phi_old = phi_new
#	x_new = f(x_old)
##				phi_new = tempphi
##				iterate+=1
##			phin[i] = phi_new
			phin[-i-1] = f(phi_old,prev,next,precision)

print phin
#		print "phi_new is ",phi_new
for i in range(len(xarray)):
	if (i!=0):
		anscheck[i] = (phin[i]-phin[i-1])*exp(phin[i]/0.026)
nconc = xarray*0
for i in range(len(xarray)):
#	if (i!=0):
	nconc[i] = exp(phin[i]/0.026)
print anscheck

#print "Local minimum occurs at ", phi_new
plt.title('Values of electron concn.(green) and phin(blue)')
plt.xlabel('length x/L')
plt.ylabel('n (in per cubic cm) and electron fermi potential (in V)')
plt.plot(xarray,phin)
#plt.plot(xarray,anscheck)
plt.plot(xarray,nconc/50)
plt.show()
