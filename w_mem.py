from numpy import *
import matplotlib.pyplot as plt
from scipy.integrate import odeint

eps = 0.01 # step size
epst = 0.001 # time step size
precision = 0.01
volt = 0.05


def f1(phi,prev,next,precision):
	value = ( prev+next*( exp(( next-prev )/0.026 )))/( 1+exp(( next-prev )/0.026 ) )
	if abs(phi-value)>precision:
		f1(value,prev,next,precision)
	return value

def f2(phi,prev,next,precision):
	value = ( prev+next*( exp(( -next+prev )/0.026 )))/( 1+exp(( -next+prev )/0.026 ) )
	if abs(phi-value)>precision:
		f2(value,prev,next,precision)
	return value


fig = plt.figure()
ax = plt.subplot(111)


volt_list = arange(0.01,volt+0.01,0.01)
#volt_list = [0.01]

for Volt in volt_list:

	EE0_e = 5.2083e10
	xarray = arange(0,1,eps)
	phin = xarray*0
	anscheck = xarray*0
	phin[-1] = Volt
	for m in range(100):
		for i in range(len(xarray)):
			if (i>0)&(i<len(xarray)-1):
				prev = phin[i-1]
				next = phin[i+1]
				phi_old = phin[i]
				phin[i] = f1(phi_old,prev,next,precision)

		for i in range(len(xarray)):
			if (i>0)&(i<len(xarray)-1):
				prev = phin[-i-2]
				next = phin[-i]
				phi_old = phin[-i-1]
				phin[-i-1] = f1(phi_old,prev,next,precision)

	for i in range(len(xarray)):
		if (i!=0):
			anscheck[i] = (phin[i]-phin[i-1])*exp(phin[i]/0.026)
	nconc = xarray*0
	for i in range(len(xarray)):
		nconc[i] = exp(phin[i]/0.026)

	nconc = nconc*(5*(10**19))

	phip = xarray*0
	anscheck = xarray*0
	phip[-1] = Volt
 

	for m in range(100):
		for i in range(len(xarray)):
			if (i>0)&(i<len(xarray)-1):
				prev = phip[i-1]
				next = phip[i+1]
				phi_old = phip[i]
				phip[i] = f2(phi_old,prev,next,precision)

		for i in range(len(xarray)):
			if (i>0)&(i<len(xarray)-1):
				prev = phip[-i-2]
				next = phip[-i]
				phi_old = phip[-i-1]
				phip[-i-1] = f2(phi_old,prev,next,precision)

	for i in range(len(xarray)):
		if (i!=0):
			anscheck[i] = (phip[i]-phip[i-1])*exp(phip[i]/0.026)
	pconc = xarray*0
	for i in range(len(xarray)):
		pconc[i] = exp(-phip[i]/0.026)

	pconc = pconc*(5*(10**19))

	ND = ones(len(xarray))
#phi = xarray*0
##ND = ND*(5*(10**19)) # ND*

	phi = phin-phip # not sure if this assumption is correct, but it holds good for the boundaries

	phi1 = phin-phip

	delnp = xarray*0

	for i in range(len(xarray)):
		if ((i>0) & (i<(len(xarray)-1))):
			delnp[i] = EE0_e * ( phi[i+1]+phi[i-1]-(2*phi[i]) )/(eps*eps) + ND[i] - (5*(10**19))


	delND = 2*precision
	delphi = 0
	iterations = 100
	NDl = ones(len(xarray))
	NDr = ones(len(xarray))
	while (delND>precision):

#ND[0] = ND[1] + (ND[1]*(phi[1]))/0.026 # From the fact that Jion[0] = 0
#ND[-1] = ND[-2] - (ND[-1]*(phi[-2]))/0.026 # From the fact that Jion[L] = 0

#for it in range(iterations):
		ND[0] = ND[1] + (ND[0]*(phi[1]))/0.026 # From the fact that Jion[0] = 0
		ND[-1] = ND[-2] - (ND[-2]*(phi[-2]))/0.026 # From the fact that Jion[L] = 0

#	while (delND>1e17):
		delND = 0

		ND_old = ND

		for i in range(len(xarray)):
			if ((i>0) & (i<(len(xarray)-1))):
				ND_oldl = ND_old[i]
		
#			ND[i] = ND[i] + 7.7e-7*epst*( (0.026*(ND[i+1]+ND[i-1]-(2*ND[i]))/(eps*eps)) + ((ND[i+1]-ND[i-1])*(phi[i+1]-phi[i-1])/(4*eps*eps)) + (ND[i]*(phi[i+1]+phi[i-1]-(2*phi[i]))/(2*eps*eps)) )
				NDl[i] = ND_old[i] + epst*( (0.026*(ND_old[i+1]+ND_old[i-1]-(2*ND_old[i]))/(eps*eps)) + ((ND_old[i+1]-ND_old[i-1])*(phi[i+1]-phi[i-1])/(4*eps*eps)) + (ND_old[i]*(phi[i+1]+phi[i-1]-(2*phi[i]))/(2*eps*eps)) )

				delND+=abs(ND_oldl - NDl[i])

		for i in range(len(xarray)):
			if ((i>0) & (i<(len(xarray)-1))):
				ND_oldr = ND_old[-i-1]
		
#			ND[-i-1] = ND[-i-1] + 7.7e-7*epst*( (0.026*(ND[-i]+ND[-i-2]-(2*ND[-i-1]))/(eps*eps)) + ((ND[-i]-ND[-i-2])*(phi[-i]-phi[-i-2])/(4*eps*eps)) + (ND[-i-1]*(phi[-i]+phi[-i-2]-(2*phi[-i-1]))/(2*eps*eps)) )
				NDr[-i-1] = ND_old[-i-1] + epst*( (0.026*(ND_old[-i]+ND_old[-i-2]-(2*ND_old[-i-1]))/(eps*eps)) + ((ND_old[-i]-ND_old[-i-2])*(phi[-i]-phi[-i-2])/(4*eps*eps)) + (ND_old[-i-1]*(phi[-i]+phi[-i-2]-(2*phi[-i-1]))/(2*eps*eps)) )

				delND+=abs(ND_oldr - NDr[-i-1])

#	print delND
		ND = (NDl+NDr)/2.0

#	print ND
#	print phi

		while (delphi>precision):
			delphi = 0
			for i in range(len(xarray)):
				if ((i>0) & (i<(len(xarray)-1))):
					phi_old = phi[i]
					phi[i] = 0.5*( phi[i+1] + phi[i-1] + ((-delnp[i]+ND[i]-(5*(10**19)))/EE0_e) )
					delphi+=abs(phi_old - phi[i])

	ND[0] = ND[1] + (ND[1]*(phi[1]))/0.026 # From the fact that Jion[0] = 0
	ND[-1] = ND[-2] - (ND[-1]*(phi[-2]))/0.026 # From the fact that Jion[L] = 0

#		delND = delND/(len(xarray)-2)

#	print phi



##def funcND(y, t):
#	return [( -y[0]/0.026 ),y[1]]
##	return (-y/0.026)

##ND0_0 = ND[0]
##ND0_1 = 1
##ND0 = [ND0_0, ND0_1]

#t = phi
#ND = odeint(funcND, ND0, phi)
##ND = odeint(funcND, ND[0], phi)
##ND_phi = xarray*0

##for i in range(len(xarray)):
##	ND_phi[i] = ND[i][0]

##print ND_phi

	delND1 = 2e17

	ND1 = ones(len(xarray))
#phi = xarray*0
#ND1 = ND1*(5*(10**19)) # ND*


	ND1[0] = (ND1[1] + (ND1[1]*(phi1[1]))/0.026) # From the fact that Jion[0] = 0
	ND1[-1] = (ND1[-1] - (ND1[-2]*(phi1[-2]))/0.026) # From the fact that Jion[L] = 0

#print "phi1 is ",phi1

#while ((delND1>1e10) & (delND1<1e25)):
#for it in range(100):
	for it in range(0):
		delND1 = 0
		for i in range(len(xarray)):
			if ((i>0) & (i<(len(xarray)-1))):
				ND_old = ND1[i]
				ND1[i] = (ND1[i-1] - ND1[i+1])*0.026*Volt/(phi1[i+1]-phi1[i-1])

				delND1+=abs(ND_old - ND1[i])

	for it in range(0):

		for i in range(len(xarray)):
			if ((i>0) & (i<(len(xarray)-1))):
				ND_old = ND1[-i-1]
				ND1[-i-1] = (ND1[-i-2] - ND1[-i])*0.026/(phi1[i+1]-phi1[i-1])

				delND1+=abs(ND_old - ND1[-i-1])



#	print delND1


	log_ND = xarray*0
	for i in range(len(xarray)):
		log_ND[i] = log(abs(ND[i]))+log(10)

#print ND


	ax.plot(xarray,log_ND,label='$v = %.2fV$'%Volt)
#	ax.plot(xarray,phi/0.026,label='$v = %.2f$'%Volt)





#for m in range(100):
#	for i in range(len(xarray)):
#		if (i!=0)&(i!=len(xarray)-1):
#			phi[i] = 0.5*( phi[i-1]+phi[i+1] + (24*eps*eps/12500)*(pconc[i]-nconc[i]+ND[i]-(5*(10**18)) ) )

ax.legend(loc=3)
plt.title('Semi-log plot of Ionic concentration - $N_{D}(x)/N_{D}^{*}$ \n (where $N_{D}^{*}=5X10^{19}cm^{-3}$)')
#plt.title('Electrostatic potential - $\phi(x)/v_{0}$ for different input voltages \n (where $v_{0}=26mV$ is the thermal voltage)')
#plt.title('Quasi-fermi potential of holes - $\phi_{p}(x)/v_{0}$ \n (where $v_{0}=%.2fV$ is input voltage)'%Volt)
#plt.title('Mobile hole concentration - $p(x)/N_{D}^{*}$ \n (where $N_{D}^{*}=5X10^{19}cm^{-3}$)')
#plt.xlabel('length x/L')
ax.set_xlabel('length x/L')
#plt.ylabel('electrostatic potential (phi) (in V)')
#ax.set_ylabel('Quasi-fermi potential of holes - $\phi_{p}(x)/v_{0}$')
#ax.set_ylabel('Mobile hole concentration - $p(x)/N_{D}^{*}$')
#ax.set_ylabel('Electrostatic potential - $\phi(x)/v_{0}$')
ax.set_ylabel('Semi-log plot of Ionic concentration - $N_{D}(x)/N_{D}^{*}$')
#plt.ylabel('(log) Ion concentration (in per cm3)')
#ax.set_ylabel('(log) Ion concentration (in per cm3)')
#plt.plot(xarray[:-2],phip[:-2])
#plt.plot(xarray,phip/Volt)
#plt.plot(xarray,anscheck)
#plt.plot(xarray,pconc/5e19)
#plt.plot(xarray,phi/0.026)
#plt.plot(xarray,ND)
#plt.plot(xarray,log_ND)
plt.show()
