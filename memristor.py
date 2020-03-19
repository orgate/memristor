from scipy.integrate import odeint
from scipy.special import gamma, airy
from numpy import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

#delt = 0.1
case = 1
timesteps = 1000
delt = 20.0/timesteps

Ron = 1
Roff = 0.5
D = 1.0
muv = 1
t = arange(0, 20.0, delt)


if (case==1):
	v = sin(t)
else:
	v = sin(t)*sin(t)
	for i in range(len(v)):
		if (i>=(len(v)/2)):
			v[i] = -1*v[i]

i = v*0
w = v*0

for n in range(len(v)):
	i[n] = v[n]/((Ron*w[n]/D) + (Roff*(1-(w[n]/D))))
	if (n<len(v)-1):
#		i[n] = v[n]/((Ron*w[n]/D) + (Roff*(1-(w[n]/D))))
		w[n+1] = w[n] + muv*(Ron/D)*i[n]*delt
#		if (w[n+1]>=1.0):
#			w[n+1] = 1.0
		#w[n] = muv*(Ron/D)*q[n]

#print v[980:]
#print i[970:]
#print w[980:]

#plt.plot(t[:960],v[:960])
#plt.plot(t[:960],i[:960])
#plt.plot(t[:980],w[:980]/D)
plt.plot(v[:960],i[:960])
plt.show()
