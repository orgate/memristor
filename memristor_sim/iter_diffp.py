from numpy import *
import matplotlib.pyplot as plt

eps = 0.01 # step size
precision = 0.00001

xarray = arange(0,1,eps)
phip = xarray*0
anscheck = xarray*0
phip[-1] = 0.04
 
def f(phi,prev,next,precision):
	value = ( prev+next*( exp(( -next+prev )/0.026 )))/( 1+exp(( -next+prev )/0.026 ) )
	if abs(phi-value)>precision:
		f(value,prev,next,precision)
	return value

for m in range(100):
	for i in range(len(xarray)):
		if (i>0)&(i<len(xarray)-1):
			prev = phip[i-1]
			next = phip[i+1]
			phi_old = phip[i]

			iterate = 0
			phip[i] = f(phi_old,prev,next,precision)

	for i in range(len(xarray)):
		if (i>0)&(i<len(xarray)-1):
			prev = phip[-i-2]
			next = phip[-i]
			phi_old = phip[-i-1]
			phip[-i-1] = f(phi_old,prev,next,precision)

print phip
for i in range(len(xarray)):
	if (i!=0):
		anscheck[i] = (phip[i]-phip[i-1])*exp(phip[i]/0.026)
nconc = xarray*0
for i in range(len(xarray)):
	nconc[i] = exp(-phip[i]/0.026)
print anscheck

plt.title('Values of hole concn.(green) and phip(blue)')
plt.xlabel('length x/L')
plt.ylabel('p (in per cubic cm) and hole fermi potential (in V)')
#plt.plot(xarray[:-2],phip[:-2])
plt.plot(xarray,phip)
#plt.plot(xarray,anscheck)
plt.plot(xarray,nconc/50)
plt.show()
