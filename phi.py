from numpy import *
import matplotlib.pyplot as plt

eps = 0.01 # step size
precision = 0.00001


EE0_e = 5.2083e10
xarray = arange(0,1,eps)
phin = xarray*0
anscheck = xarray*0
phin[-1] = 0.05

def f(phi,prev,next,precision):
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
			phin[i] = f(phi_old,prev,next,precision)

	for i in range(len(xarray)):
		if (i>0)&(i<len(xarray)-1):
			prev = phin[-i-2]
			next = phin[-i]
			phi_old = phin[-i-1]
			phin[-i-1] = f(phi_old,prev,next,precision)

for i in range(len(xarray)):
	if (i!=0):
		anscheck[i] = (phin[i]-phin[i-1])*exp(phin[i]/0.026)
nconc = xarray*0
for i in range(len(xarray)):
	nconc[i] = exp(phin[i]/0.026)

nconc = nconc*(2*(10**19))

phip = xarray*0
anscheck = xarray*0
phip[-1] = 0.05
 
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
			phip[i] = f(phi_old,prev,next,precision)

	for i in range(len(xarray)):
		if (i>0)&(i<len(xarray)-1):
			prev = phip[-i-2]
			next = phip[-i]
			phi_old = phip[-i-1]
			phip[-i-1] = f(phi_old,prev,next,precision)

for i in range(len(xarray)):
	if (i!=0):
		anscheck[i] = (phip[i]-phip[i-1])*exp(phip[i]/0.026)
pconc = xarray*0
for i in range(len(xarray)):
	pconc[i] = exp(-phip[i]/0.026)

pconc = pconc*(2*(10**19))

ND = ones(len(xarray))
#phi = xarray*0
ND = ND*(5*(10**19)) # ND*

phi = phin-phip # not sure if this assumption is correct, but it holds good for the boundaries

#for m in range(100):
#	for i in range(len(xarray)):
#		if (i!=0)&(i!=len(xarray)-1):
#			phi[i] = 0.5*( phi[i-1]+phi[i+1] + (24*eps*eps/12500)*(pconc[i]-nconc[i]+ND[i]-(5*(10**18)) ) )

plt.title('Values of phi')
plt.xlabel('length x/L')
plt.ylabel('electrostatic potential (phi) (in V)')
#plt.plot(xarray[:-2],phip[:-2])
#plt.plot(xarray,phip)
#plt.plot(xarray,anscheck)
#plt.plot(xarray,pconc/50)
plt.plot(xarray,phi)
plt.show()
