import ROOT
import numpy as np
from matplotlib import pyplot as plt
import array
import os, sys
from multiprocessing import Pool

rate, efficiency = np.loadtxt("MVA-points.txt", unpack = True)

# DATA GENERATION
simulationCPP = "LoopLineShape"
mvaScan = "\"true\""

ConfFile = "\"ToyConfiguration.txt\""


if "-g" in sys.argv:
	ROOT.gInterpreter.ProcessLine(".L LoopLineShape.cpp")
	Nfiles = sys.argv[2]
	for i, item in enumerate(efficiency):
		CosmicRate = str(rate[i])
		Efficiency = str(item)
		folder = "mva_" + str(i)
		# Create folder
		
		if not os.path.exists(folder):	# Check if the folder exist
			os.makedirs(folder)	# Create the folder

		code = simulationCPP + "(" + Nfiles + ",\"" + folder + "/\"," + mvaScan + "," + CosmicRate + "," + Efficiency + "," + ConfFile + ")"
		print("running: ", code)
		# Execute simulation code
		ROOT.gInterpreter.ProcessLine(code)
		# copy the configuration file in the folder
		os.popen("cp ToyConfiguration.txt " + folder + "/"  + ConfFile)


# Constant Fraction
bias_cf = np.array([])
variance_cf = np.array([])
# Foward 2017
bias_fw = np.array([])
variance_fw = np.array([])
# Reversed 2017
bias_rev = np.array([])
variance_rev = np.array([])
# Sum Neighbors 
bias_sum = np.array([])
variance_sum = np.array([])
# Threshold
bias_thr = np.array([])
variance_thr = np.array([])

# DATA ANALYSIS
# Execute analysis code


ROOT.gInterpreter.ProcessLine(".L ScanAnalysis.cpp")
folder = "mva_"
trial = 199

def task(i):
	values = ROOT.ScanAnalysis(folder + str(i) + "/", folder +  str(i) + "/ToyConfiguration.txt",0,trial,3,0.3,rate[i])
	npResult = np.asarray(values)
	return npResult

uf __name__ = "__main__":
	# create the process pool
	index = range(0,len(rate))
	with Pool(processes = 8) as pool:
		for results in pool.map(task, index)
			np.append(bias_thr, results[0])
			np.append(variance_thr, results[1])
			np.append(bias_fw, results[2])
			np.append(variance_fw, results[3])
			np.append(bias_rev, results[4])
			np.append(variance_rev, result[5])
			np.append(bias_cf, results[6])
			np.append(variance_cf, results[7])
			np.append(bias_sum, results[8])
			np.append(variance_sum, results[9])

"""
#SEQUENTIAL VERSION

# Constant Fraction
bias_cf = np.array(len(rate))
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

for i, item in enumerate(rate):
	result = ROOT.ScanAnalysis(folder + str(i) + "/", folder +  str(i) + "/ToyConfiguration.txt",0,trial,3,0.3,item)
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
	bias_sum[i] = npResult[8]
	variance_sum[i] = npResult[9]
"""

cost_cf = (np.abs(bias_cf) + variance_cf)**2/variance_cf	
cost_fw = (np.abs(bias_fw) + variance_fw)**2/variance_fw
cost_rev = (np.abs(bias_rev) + variance_rev)**2/variance_rev
cost_sum = (np.abs(bias_sum) + variance_sum)**2/variance_sum
cost_thr = (np.abs(bias_thr) + variance_thr)**2/variance_thr

plt.figure(1)
plt.title("MVA scan " + r"N_{trial} = " + str(trial + 1) , fontsize = 20)
plt.xlabel("rate [events/second]", fontsize = 15)
plt.ylabel("cost [kHz]", fontsize = 15)
plt.xscale("log")
plt.grid()
plt.errorbar(rate,cost_cf, linestyle = '',marker = "v", color = "red", label = "constant fraction 30%")
plt.errorbar(rate,cost_fw, linestyle = '',marker = ".", color = "blue", label = "foward")
plt.errorbar(rate,cost_sum, linestyle = '',marker = "s", color = "green", label = "Neighbors sum (N = 3)")
plt.errorbar(rate,cost_thr, linestyle = '',marker = "D", color = "brown", label = "over Threshold (> 3)")
plt.legend(fontsize = 10)

plt.figure(2)
plt.title("Bias " + r"N_{trial} = " + str(trial + 1) , fontsize = 20)
plt.xlabel("rate [events/second]", fontsize = 15)
plt.ylabel("Bias", fontsize = 15)
plt.xscale("log")
plt.grid()
plt.errorbar(rate,bias_cf, linestyle = '',marker = "v", color = "red", label = "constant fraction 30%")
plt.errorbar(rate,bias_fw, linestyle = '',marker = ".", color = "blue", label = "foward")
plt.errorbar(rate,bias_sum, linestyle = '',marker = "s", color = "green", label = "Neighbors sum (N = 3)")
plt.errorbar(rate,bias_thr, linestyle = '',marker = "D", color = "brown", label = "over Threshold (> 3)")
plt.legend(fontsize = 10)

plt.figure(3)
plt.title("Variance vs Bias " + r"N_{trial} = " + str(trial + 1) , fontsize = 20)
plt.xlabel("Variance [kHz]", fontsize = 15)
plt.ylabel("Bias [kHz]", fontsize = 15)
plt.grid()
plt.errorbar(variance_cf,bias_cf, linestyle = '',marker = "v", color = "red", label = "constant fraction 30%")
plt.errorbar(variance_fw,bias_fw, linestyle = '',marker = ".", color = "blue", label = "foward")
plt.errorbar(variance_sum,bias_sum, linestyle = '',marker = "s", color = "green", label = "Neighbors sum (N = 3)")
plt.errorbar(variance_thr,bias_thr, linestyle = '',marker = "D", color = "brown", label = "over Threshold (> 3)")
plt.legend(fontsize = 10)

plt.show()
