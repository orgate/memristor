# From calculation, we expect that the local minimum occurs at x=9/4

from numpy import *
import matplotlib.pyplot as plt

#x_old = 0
#x_new = 2 # The algorithm starts at x=6
eps = 0.01 # step size
precision = 0.01

xarray = arange(0,1,eps)
phin = xarray*0
anscheck = xarray*0
phin[-1] = 0.1
length = len(xarray)
#def f_prime(x):
#    return 4 * x**3 - 9 * x**2

##def f(phi,prev,next):
##	value = ( prev+next*( exp(( next-phi )/0.026 )))/( 1+exp(( next-phi )/0.026 ) )
#	value = -exp(x)
##	return value
 
def f(phi,prev,next,precision,i):
#	value = ( prev+next*( exp(( next-phi )/0.026 )))/( 1+exp(( next-phi )/0.026 ) )
	value = phi - ( (((next[0]-prev[1])/(2*eps))**2)/0.026 + (next[0]+prev[1]-(2*phi))/(eps*eps) )/( (((next[0]-prev[1])/(2*eps))**2)/(0.026**2) + (3/0.026)*(next[0]+prev[1]-(2*phi))/(eps*eps) + ( (prev[0]- 2*prev[1] + 2*next[0] - next[1])/(2*(eps**3)) )/( (next[0]-prev[1])/(2*eps) ) )
	if abs(phi-value)>precision:
		f(value,prev,next,precision,i)
	return value

for m in range(20):
	for i in range(len(xarray)):
		if (i>0)&(i<len(xarray)-1):
			if (i==1):
				prev = [phin[-i-3],phin[-i-2]]
				next = [phin[-i],phin[-1]]
			elif (i==length-2):
				prev = [phin[0],phin[-i-2]]
				next = [phin[-i],phin[-i+1]]
			else:
				prev = [phin[-i-3],phin[-i-2]]
				next = [phin[-i],phin[-i+1]]
			phi_old = phin[-i-1]
##			phi_new = 10	

##			while (abs(phi_new - phi_old) > precision)&(iterate<10):
##				tempphi = f(phi_old,prev,next)
##				phi_old = phi_new
#	x_new = f(x_old)
##				phi_new = tempphi
##				iterate+=1
##			phin[i] = phi_new
			phin[-i] = f(phi_old,prev,next,precision,(length-1-i))

	for i in range(len(xarray)):
		if (i>0)&(i<len(xarray)-1):
			if (i==1):
				prev = [phin[0],phin[i-1]]
				next = [phin[i+1],phin[i+1]]
			elif (i==length-2):
				prev = [phin[i-2],phin[i-1]]
				next = [phin[i+1],phin[-1]]
			else:
				prev = [phin[i-2],phin[i-1]]
				next = [phin[i+1],phin[i+2]]
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
			phin[i] = f(phi_old,prev,next,precision,i)

#print phin


print phin
#		print "phi_new is ",phi_new
for i in range(len(xarray)):
	if (i!=0):
		anscheck[i] = (phin[i]-phin[i-1])*exp(phin[i]/0.026)
print anscheck

#print "Local minimum occurs at ", phi_new
plt.plot(xarray,phin)
#plt.plot(xarray,anscheck)
plt.show()
