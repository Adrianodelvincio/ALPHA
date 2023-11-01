import ROOT
from ROOT import *
from ROOT import TMath
from array import array

#defining the fit function

def myRayleigh(x, p):
    sigma = par[0]
    r = x
    tmp = r**2/sigma * np.exp(-(x)**2/(2*sigma))
    return tmp

#f1 = ROOT.TF1("funzione", myRayleigh, 0,4)
f1 = ROOT.TF1("funzione","(x/[0]**2)*exp((-(x**2/(2*[0]**2))))" ,0.0001,4)

par = array('d', [1])
print(par)
f1.SetParameters(par)
c = ROOT.TCanvas()
c.Draw() # Necessary to make the graphics show!
f1.Draw("SAME")
input()
