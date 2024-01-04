import ROOT
import numpy as np

result = ROOT.gInterpreter.ProcessLine(".x code.cpp")
npresult = np.asarray(result)
print( type(npresult), "value: ", npresult)
