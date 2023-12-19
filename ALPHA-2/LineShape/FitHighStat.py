import ROOT
import numpy as np

window,freq,power,unixtimestart,runtimestart,stop,duration,Trig,Read,Pass,MVA = np.loadtxt("ToyModelData.txt", unpack = True)

freq = freq * 1e6 - 28.2353e6 	# Eliminate the offset and convert to kHz
background = 0.051028571*8		# Expected background rates
print("max frequence: ", np.max(freq), "number of bin: " ,len(freq))

h = ROOT.TH1F('lineshape', 'highStat', len(freq), -2.5 , np.max(freq) - 2.5)

for i, item in enumerate(Pass):
	h.SetBinContent(i+1,item)

baseline = ROOT.TF1('f1', '[0]', 0, 150)	# define the fit function, a baseline
baseline.SetRange(0, 100)					# set the range of the function
h.Fit(baseline, 'R')						# fit the histogram using the range of the function
baseline.SetParNames('baseline')
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)

c1 = ROOT.TCanvas("c1", "c1",800,800)
h.SetAxisRange(0,40, "Y")
#h.SetAxisRange(-40,320, "X")
h.SetLineWidth(3)
h.SetMarkerStyle(23);
h.SetMarkerColor(1)
h.SetMarkerSize(1)
h.DrawClone("")
h.DrawClone("SAME P")
c1.Print("Fit.pdf")

mask = (freq <= 50)
print("mean is: ", np.mean(Pass[mask]))
