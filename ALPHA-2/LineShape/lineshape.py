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
	if(f < x0):
		arg = np.exp(-(f-x0)**2/(2*sigma0**2 + k0*(f-x0)**2))
	if(f > x0):
		arg = np.exp(-(f-x0)**2/(2*sigma1**2 + k1*(f-x0)**2))
	return N*arg

Cruijff = np.vectorize(Cruijff)
peak = 270 ; xpeak = 220
onset = 175

plt.figure(1)
# bellurie
plt.grid()
plt.title("LineShape")
plt.xlabel("kHz")
plt.ylabel("Counts")
#plot the data
plt.errorbar(freq,Pass, linestyle = '', marker = 's', color = 'red', markersize = 3)
plt.step(freq,Pass, linestyle = '-', color = 'black', where='mid')
#plot the functions
xx = np.linspace(150, xpeak,100)
plt.plot(xx, linearRise(xx, onset, peak/(xpeak - onset) ), linestyle = '--', color = 'purple', label = "linear rise")
#xx = np.linspace(onset, xpeak, 100)
xx = np.linspace(175, xpeak,100)
plt.plot(xx, quadratic(xx, onset, peak/(xpeak- onset)**2), linestyle = '--', color = 'green', label = "quadratic")

onset = 150 ; xpeak = 215 ; peak = 230

plt.plot(xx, expRise(xx, onset, np.log(peak)/(xpeak - onset)))
plt.legend()


plt.figure(2)
plt.errorbar(freq,Pass, linestyle = '', marker = 's', color = 'red', markersize = 3)
plt.step(freq,Pass, linestyle = '-', color = 'black', where='mid')
Nnorm = 110#Nnorm = np.sum(Pass);
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

plt.figure(3)
plt.xlabel("frequency [kHz]")
plt.ylabel("counts")
plt.xlim(100,320)
plt.errorbar(freq,Pass, linestyle = '', marker = 's', color = 'black', markersize = 3)
plt.step(freq,Pass, linestyle = '-', color = 'black', where='mid')
xx = np.linspace(100, 600,1000)
plt.plot(xx,Cruijff(xx,225,12,38,0.1,0,270), linestyle = '--', color = 'red', label = r'alpha = %.1f ; beta = %.1f' % (alpha, b))
plt.show()
