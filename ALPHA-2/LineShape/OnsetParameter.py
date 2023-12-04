import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit as fit
import scipy.integrate as integrate

def Cruijff(f,x0,sigma0, sigma1, k0, k1,N):
	arg = 0
	if(f < x0):
		arg = np.exp(-(f-x0)**2/(2*sigma0**2 + k0*(f-x0)**2))
	if(f >= x0):
		arg = np.exp(-(f-x0)**2/(2*sigma1**2 + k1*(f-x0)**2))
	return N*arg

Cruijff = np.vectorize(Cruijff)

sigma0 = 8.8
sigma1 = 40.7
x0 = 220
k0 = 0.22
k1 = -0.03
N = 271

result = integrate.quad(lambda x: Cruijff(x,x0,sigma0, sigma1, k0, k1,N), 0, 316)
print("result of the integration: ", result[0] , " +/- ", result[1])

plt.figure(1)
plt.title("Cruijff")
plt.grid()
start = 0
stop = 320
xx = np.linspace(start, stop, 1000)
plt.plot(xx,(1/result[0])*Cruijff(xx,x0,sigma0,sigma1,k0,k1,N), linestyle = '--', color = 'blue', label = "Cruijff")

check = integrate.quad(lambda x: (1/result[0])*Cruijff(x,x0,sigma0, sigma1, k0, k1,N), 0, 316)
print("let's see if it is ok: ", check[0])
Nstack = 20; Hbar = 14
Ntot = Nstack*Hbar
Nstep = 19
print("total number anti-hydrogen: ", Ntot)


x = np.linspace(start,stop,Nstep) # Defining the swipe
# Binwidth * pdf at center normalized
y = (np.diff(x)[0]) * (1/result[0])*Cruijff(x,x0,sigma0,sigma1,k0,k1,N)
y = y/y.sum() 	# Normalizzo la lineshape campionata
y = Ntot*y 	# Scalo con la statistica in esame
print("Ma y è normalizzato giusto?: ", y.sum(), " Ntot: ", Ntot)

leftIntegral = np.empty(len(y))
somma = 0
for i,item in enumerate(y):
	somma += item 
	leftIntegral[i] = somma

mask = (y > 0.9)
mask2 = (leftIntegral > 1)
plt.figure(2, figsize = (15,10))
plt.grid()
plt.title("LineShape for run 4b")
plt.plot(x,y, linestyle = '-', marker = '.' ,color = 'black', label = "sampling")
plt.axvline(x = min(x[mask]), ymin=0, ymax=10, label = "onset at N expected = 1", linestyle = '--' ,color = 'red')
plt.axvline(x = min(x[mask2]), ymin=0, ymax=10, label = r"onset at N($f < \overline{f}$ = 1)", linestyle = '--' ,color = 'purple')
plt.legend()
plt.show()
