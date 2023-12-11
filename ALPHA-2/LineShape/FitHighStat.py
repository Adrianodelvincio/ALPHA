import ROOT
import numpy as np

window,freq,power,unixtimestart,runtimestart,stop,duration,Trig,Read,Pass,MVA = np.loadtxt("ToyModelData.txt", unpack = True)

freq = freq * 1e6 - 28.2353e6 # Eliminate the offset and convert to kHz
background = 0.051028571*8
print("max frequence: ", np.max(freq), "number of bin: " ,len(freq))

h = ROOT.TH1F('lineshape', 'highStat', len(freq), -2.5 , np.max(freq) - 2.5)

for i, item in enumerate(Pass):
	h.SetBinContent(i,item)


c1 = ROOT.TCanvas("c1", "c1",800,800)
h.DrawClone()
c1.Print("Fit.pdf")
