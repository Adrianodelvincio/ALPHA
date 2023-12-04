import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit as fit
from scipy.stats import beta

window,freq,power,unixtimestart,runtimestart,stop,duration,Trig,Read,Pass,MVA = np.loadtxt("elog.txt", unpack = True)

freq = freq * 1e6 - 28.2353e6

def linearRise(f,f0,m):
	tmp = (f-f0)*m
	if(f < f0):
		tmp = 0
	return tmp

linearRise = np.vectorize(linearRise)

def expRise(f,f0,m):
	return np.exp(m*(f - f0))

def quadratic(f, f0, m):
	return m*(f-f0)**2
	
def Cruijff(f,x0,sigma0, sigma1, k0, k1,N):
	arg = 0
	if(f < x0):
		arg = np.exp(-(f-x0)**2/(2*sigma0**2 + k0*(f-x0)**2))
	if(f >= x0):
		arg = np.exp(-(f-x0)**2/(2*sigma1**2 + k1*(f-x0)**2))
	return N*arg

Cruijff = np.vectorize(Cruijff)


plt.figure(1)
# bellurie
plt.grid()
plt.title("LineShape")
plt.xlabel("kHz")
plt.ylabel("Counts")
#plot the data
plt.errorbar(freq,Pass, linestyle = '', marker = 's', color = 'red', markersize = 3)
plt.step(freq,Pass, linestyle = '-', color = 'black', where='mid')

###plot the functions
peak = 270 ; xpeak = 220
onset = 175 ; xx = np.linspace(150, xpeak,100)
#linear
plt.plot(xx, linearRise(xx, onset, peak/(xpeak - onset)), linestyle = '--', color = 'purple', label = "linear rise")
xx = np.linspace(onset, xpeak,100)
#quadratic
plt.plot(xx, quadratic(xx, onset, peak/(xpeak- onset)**2), linestyle = '--', color = 'green', label = "quadratic")
#exponential
onset = 150 ; xpeak = 215 ; peak = 230
plt.plot(xx, expRise(xx, onset, np.log(peak)/(xpeak - onset)))
plt.legend()
# Beta distribution
plt.figure(2)
plt.errorbar(freq,Pass, linestyle = '', marker = 's', color = 'red', markersize = 3)
plt.step(freq,Pass, linestyle = '-', color = 'black', where='mid')
Nnorm = 110	#Nnorm = np.sum(Pass);
onset = 175
xx = np.linspace(onset, 316, 100)

alpha = 2; b = 5;
plt.plot(xx, Nnorm*beta.pdf((xx - onset)/(xx.max() - onset), alpha, b), linestyle = '--', label = r'alpha = %.1f ; beta = %.1f' % (alpha, b))
alpha = 2.5; b = 5;
plt.plot(xx, Nnorm*beta.pdf((xx - onset)/(xx.max() - onset), alpha, b), linestyle = '--', label = r'alpha = %.1f ; beta = %.1f' % (alpha, b))
alpha = 3; b = 5;
plt.plot(xx, Nnorm*beta.pdf((xx - onset)/(xx.max() - onset), alpha, b), linestyle = '--', label = r'alpha = %.1f ; beta = %.1f' % (alpha, b))
alpha = 4; b = 5;
plt.plot(xx, Nnorm*beta.pdf((xx - onset)/(xx.max() - onset), alpha, b), linestyle = '--', label = r'alpha = %.1f ; beta = %.1f' % (alpha, b))
alpha = 4; b = 6;
plt.plot(xx, Nnorm*beta.pdf((xx - onset)/(xx.max() - onset), alpha, b), linestyle = '--', label = r'alpha = %.1f ; beta = %.1f' % (alpha, b))
alpha = 4; b = 7;
plt.plot(xx, Nnorm*beta.pdf((xx - onset)/(xx.max() - onset), alpha, b), linestyle = '--', label = r'alpha = %.1f ; beta = %.1f' % (alpha, b))
alpha = 4; b = 4;
plt.plot(xx, Nnorm*beta.pdf((xx - onset)/(xx.max() - onset), alpha, b), linestyle = '--', label = r'alpha = %.1f ; beta = %.1f' % (alpha, b))
alpha = 4; b = 8;
plt.plot(xx, Nnorm*beta.pdf((xx - onset)/(xx.max() - onset), alpha, b), linestyle = '--', label = r'alpha = %.1f ; beta = %.1f' % (alpha, b))
plt.legend()

#cruijff function
x0 = 225
sigma0 = 12
sigma1 = 38
k0 = 0.1
k1 = 0
N = 270
fig, (ax1, ax2) = plt.subplots(2, figsize = (15,9), height_ratios=[2, 1] , layout = 'tight')
ax1.grid()
ax1.set_title("Cruijff fit to LineShape", fontsize = 18, color = 'blue')
#ax1.set_xlabel("frequency [kHz]")
ax1.set_ylabel("Counts", fontsize = '12')
ax1.set_xlim(100,320)
ax1.errorbar(freq,Pass, linestyle = '', marker = 's', color = 'black', markersize = 3)
ax1.step(freq,Pass, linestyle = '-', color = 'black', where='mid')
xx = np.linspace(100, 320,1000)
#ax1.plot(xx,Cruijff(xx,x0,sigma0,sigma1,k0,k1,N), linestyle = '--', color = 'red', label = r"$\sigma_{0} = %.1f$" "\n" r"$\sigma_{1} = %.1f$" "\n" r"$x_{0} = %.1f$" % (sigma0, sigma1, x0))
## Fit to the data with Cruijff function
mask = (Pass >= 3)
popt, pcovm = fit(Cruijff, freq[mask], Pass[mask], p0 = [x0,sigma0,sigma1,k0,k1,N])

mask = (Pass > 3)
chisq = ((Pass[mask] - Cruijff(freq[mask],*popt))**2/(Pass[mask])).sum()

ax1.plot(xx,Cruijff(xx,*popt), linestyle = '--', color = 'red', label = r"$\sigma_{0} = %.1f$" "\n" r"$\sigma_{1} = %.1f$" "\n" r"$x_{0} = %.1f$" "\n" r"$k_{0} = %.2f$" "\n" r"$k_{1} = %.2f$" "\n" r"$N = %.1f$" "\n" r"$ \frac{\chi^{2}}{ndof} = \frac{%.1f}{%d} \pm %.1f$" % (popt[1], popt[2], popt[0], popt[3], popt[4], popt[5], chisq, len(Pass[mask]), np.sqrt(2*len(Pass[mask]))))
ax1.legend( fontsize = '13')

# residual
ax2.grid()
ax2.set_xlim(100,320)
ax2.set_title("Residuals", fontsize = '12')
ax2.set_xlabel("frequency [kHz]", fontsize = '12')
ax2.set_ylabel("Counts", fontsize = '12')
Residui = Pass - Cruijff(freq,x0,sigma0,sigma1,k0,k1,N)
ax2.errorbar(freq, Residui,marker = '.', linestyle = 'dotted', color = 'green')
fig.savefig("Plot/FitToLineShape.pdf", format = 'pdf' , bbox_inches = 'tight')
plt.show()

