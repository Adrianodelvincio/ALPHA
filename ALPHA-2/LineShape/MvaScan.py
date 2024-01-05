import ROOT
import numpy as np
from matplotlib import pyplot as plt
import array
import os

rate, efficiency = np.loadtxt("MVA-points.txt", unpack = True)

# DATA GENERATION
simulationCPP = "LoopLineShape"
mvaScan = "\"true\""
Nfiles = "30"
ConfFile = "\"ToyConfiguration.txt\""


ROOT.gInterpreter.ProcessLine(".L LoopLineShape.cpp")
for i, item in enumerate(efficiency):
	CosmicRate = str(rate[i])
	Efficiency = str(item)
	folder = "mva_" + str(i)
	# Create folder	
	if not os.path.exists("/home/commodo98/Documenti/ALPHA/ALPHA-2/LineShape/" + folder):	# Check if the folder exist
		os.makedirs("/home/commodo98/Documenti/ALPHA/ALPHA-2/LineShape/" + folder)			# Create the folder		
	
	code = simulationCPP + "(" + Nfiles + ",\"" + folder + "/\"," + mvaScan + "," + CosmicRate + "," + Efficiency + "," + ConfFile + ")"
	print("running: ", code)
	# Execute simulation code
	ROOT.gInterpreter.ProcessLine(code)
	# copy the configuration file in the folder
	os.popen("cp ToyConfiguration.txt " + folder + "/"  + ConfFile)


# Constant Fraction
bias_cf = np.zeros(len(rate))
variance_cf = np.zeros(len(rate))
# Foward 2017
bias_fw = np.zeros(len(rate))
variance_fw = np.zeros(len(rate))
# Reversed 2017
bias_rev = np.zeros(len(rate))
variance_rev = np.zeros(len(rate))
# Sum Neighbors 
bias_sum = np.zeros(len(rate))
variance_sum = np.zeros(len(rate))
# Threshold
bias_thr = np.zeros(len(rate))
variance_thr = np.zeros(len(rate))

# DATA ANALYSIS
# Execute analysis code
ROOT.gInterpreter.ProcessLine(".L ScanAnalysis.cpp")
for i, item in enumerate(rate):
	result = ROOT.std.vector("double")(10)
	result = ROOT.ScanAnalysis("mva_" + str(i) + "/", "mva_" +  str(i) + "/ToyConfiguration.txt",0,30,3,0.1,item)
	print("point %d" % i, " ", type(result), " ", result)
	npResult = np.asarray(result)
	
	bias_thr[i] = npResult[0]
	variance_thr[i] = npResult[1]
	bias_cf[i] = npResult[6]
	variance_cf[i] = npResult[7]
	bias_fw[i] = npResult[2]
	variance_fw[i] = npResult[3]
	bias_rev[i] = npResult[4]
	variance_rev[i] = npResult[5]
	bias_sum = npResult[8]
	variance_sum[i] = npResult[9]
	
plt.figure(1)
plt.title("MVA scan", fontsize = 20)
plt.xlabel("rate [events/second]")
plt.ylabel("Variace")
plt.xscale("log")
plt.grid()
plt.errorbar(rate,variance_cf, linestyle = '',marker = "v", color = "red", label = "constant fraction 10%")
plt.errorbar(rate,variance_fw, linestyle = '',marker = ".", color = "blue", label = "foward")
plt.errorbar(rate,variance_sum, linestyle = '',marker = "s", color = "green", label = "Neighbors sum")
plt.errorbar(rate,variance_thr, linestyle = '',marker = "D", color = "brown", label = "over Threshold (> 3)")
plt.legend()
plt.show()
