import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit as fit
import scipy.integrate as integrate

def Cruijff(f,onset,x0,sigma0, sigma1, k0, k1,N):
	arg = 0
	if  (f < onset):
		arg = 0;
	elif(f < x0):
		arg = np.exp(-(f-x0)**2/(2*sigma0**2 + k0*(f-x0)**2))
	elif(f >= x0):
		arg = np.exp(-(f-x0)**2/(2*sigma1**2 + k1*(f-x0)**2))
	return N*arg

Cruijff = np.vectorize(Cruijff)

# Parameters of the Cruijff
x_cb_start = 175
x0 = 220.58
sigma0 = 9.06
sigma1 = 40.49
k0 = 0.20
k1 = -0.03
Norm = 271.16

## CHECK THE NORMALIZATION OF THE SAMPLEDLINESHAPE
Nstack = 200; Hbar = 14
Ntot = Nstack*Hbar/2
Nstep = 50
SweepStep = 24

start = x_cb_start - 3*5
stop = start + Nstep*5

x = np.arange(start,stop,5) # Total Number of Swipes
y = Cruijff(x,x_cb_start,x0,sigma0,sigma1,k0,k1,Norm)
y = y/y.sum() 	# Normalizzo la lineshape campionata
print("Normalized peak values of the lineshape: ", np.max(y))
y = Ntot*y 	# Scalo con la statistica in esame
print("Total Number of event sampled: ", int(y.sum()), " Ntot: ", Ntot)

mask = (x < start + 5*SweepStep)
Hbar_effective = (y[mask]).sum()
print("Considerando 24 SweepStep, numero di eventi: ", Hbar_effective)
print("Percentuale: ", Hbar_effective/Ntot)

plt.figure(0)
plt.grid()
plt.errorbar(x,y, marker = 's', linestyle = '', color = 'red', label = 'total sampled')
plt.errorbar(x[mask], y[mask], linestyle = 'dotted', color = 'blue')
plt.legend()
plt.show()



### Now do the fit with a step-function

window,freq,power,unixtimestart,runtimestart,stop,duration,Trig,Read,Pass,MVA = np.loadtxt("ToyModelData.txt", unpack = True)

freq = freq * 1e6 - 28.2353e6 # Eliminate the offset and convert to kHz
background = 0.051028571*8
#Pass = Pass - background

mask0 = (Pass >= 1)
mask1 = (freq >= 100)
mask = mask0 & mask1

def Cruijff(f,x0,sigma0, sigma1, k0, k1,N,baseline):
	arg = 0
	if(f < x0):
		arg = N*np.exp(-(f-x0)**2/(2*sigma0**2 + k0*(f-x0)**2))
	if(f >= x0):
		arg = N*np.exp(-(f-x0)**2/(2*sigma1**2 + k1*(f-x0)**2))
	onset = 175
	if( f <= onset):
		arg = baseline
	return arg

Cruijff = np.vectorize(Cruijff)

yerr = np.sqrt(Pass)
mask2 = (Pass == 0); mask3 = mask2 & mask1
 
yerr[mask3] = 1
popt, pcovm = fit(Cruijff, freq[mask1], Pass[mask1], p0 = [x0,sigma0,sigma1,k0,k1,Norm,1], sigma = yerr[mask1])

#Print the result of the fit
print("parametri: ", popt)
print("errori:    ", np.sqrt(pcovm.diagonal()))
errors = np.sqrt(pcovm.diagonal())
chisq = ((Pass[mask] - Cruijff(freq[mask],*popt))**2/(Pass[mask])).sum()

#save the result on the Configuration file
configFile = open("ToyConfiguration.txt", "a")
configFile.write("\n" "x0 = %.2f" % popt[0])
configFile.write("\n" "sigma0 = %.2f" % popt[1])
configFile.write("\n" "sigma1 = %.2f" % popt[2])
configFile.write("\n" "k0 = %.2f" % popt[3])
configFile.write("\n" "k1 = %.2f" % popt[4])
configFile.write("\n" "Norm = %.2f" % popt[5])

configFile.close()

#Plot the lineshape and the cruijff
fig, (ax1, ax2) = plt.subplots(2, figsize = (15,9), height_ratios=[2, 1] , layout = 'tight')
ax1.grid()
ax1.set_title("Cruijff fit to LineShape", fontsize = 18, color = 'blue')
ax1.set_ylabel("Counts", fontsize = '12')
ax1.set_xlim(100,320)
ax1.errorbar(freq,Pass,yerr, capsize = 3  ,linestyle = '', marker = 's', color = 'blue', markersize = 3)
ax1.step(freq,Pass, linestyle = '-', color = 'black', where='mid')
xx = np.linspace(100, 320,1000)
ax1.plot(xx,Cruijff(xx,*popt), linestyle = '--', color = 'red', label = r"$\sigma_{0} = %.1f \pm %.2f$" "\n" r"$\sigma_{1} = %.1f \pm %.2f$" "\n" r"$x_{0} = %.1f \pm %.2f$" "\n" r"$k_{0} = %.2f \pm %.2f$" "\n" r"$k_{1} = %.2f \pm %.2f$" "\n" r"$N = %.1f \pm %.2f$" "\n" r"baseline $= %.1f \pm %.1f$" "\n" r"$ \frac{\chi^{2}}{ndof} = \frac{%.1f}{%d} \pm %.1f$" % (popt[1], errors[1] ,popt[2], errors[2], popt[0], errors[0], popt[3], errors[3], popt[4], errors[4], popt[5],errors[5],popt[6], errors[6], chisq, len(Pass[mask]), np.sqrt(2*len(Pass[mask]))))
ax1.legend( fontsize = '13')
# residual
ax2.grid()
ax2.set_xlim(100,320)
ax2.set_title("Residuals", fontsize = '12')
ax2.set_xlabel("frequency [kHz]", fontsize = '12')
ax2.set_ylabel("Normalized Residual", fontsize = '12')
Residui = (Pass[mask1] - Cruijff(freq[mask1],*popt))/(yerr[mask1])
ax2.errorbar(freq[mask1], Residui,marker = '.', linestyle = 'dotted', color = 'green')
fig.savefig("Plot/TruncatedLineShape.pdf", format = 'pdf' , bbox_inches = 'tight')
plt.show()

