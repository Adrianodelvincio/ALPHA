import ROOT
import numpy as np
from matplotlib import pyplot as plt
import array
import os, sys
from multiprocessing import Pool
from matplotlib import rcParams
rcParams['axes.titlepad'] = 20 

rate, efficiency = np.loadtxt("MVA-points.txt", unpack = True)

# DATA GENERATION
simulationCPP = "LoopLineShape"
mvaScan = "\"true\""
ConfFile = "\"ToyConfiguration.txt\""

def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

if "-g" in sys.argv:
	Nfiles = sys.argv[2]
	ROOT.gInterpreter.ProcessLine(".L LoopLineShape.cpp")

	def generate(i):
		CosmicRate = str(rate[i])
		Efficiency = str(efficiency[i])
		folder = "try_" + str(i)
		if not os.path.exists(folder):	# Check if the folder exist
			os.makedirs(folder)	# Create the folder
		code = simulationCPP + "(" + Nfiles + ",\"" + folder + "/\"," + mvaScan + "," + CosmicRate + "," + Efficiency + "," + ConfFile + ")"
		print("running: ", code)
		ROOT.LoopLineShape(int(Nfiles),folder + "/","true", rate[i], efficiency[i],"ToyConfiguration.txt")
		os.popen("cp ToyConfiguration.txt " + folder + "/"  + "ToyConfiguration.txt")	
	
	# Paraller generation of the data
	njob = 7
	indexes = range(0,len(rate))
	for i in range(0, int(len(rate)/njob)):
		#Create list of point to be generated
		start = i*njob
		stop = (i+1)*njob
		index = np.asarray(indexes[start:stop])
		print("Generating MVA points: ", index)
		if __name__ == "__main__":
			# Generate Files
			with Pool(processes = njob, maxtasksperchild = 1) as pool:
				pool.map(generate, index)
		# Create list of last points to be generated
		if (i == (int(len(rate)/njob) - 1) and (len(rate))%njob != 0):
			start = (i+1)*njob
			stop = start + len(rate)%njob
			index = np.asarray(indexes[start:stop])
			print("generating last points: " , index)
			# Generate last points
			if __name__ == "__main__":
				with Pool(processes = njob, maxtasksperchild = 1) as pool:
					pool.map(generate, index)
	for i in range(0,len(rate)):
		folder = "mva_" + str(i)
		replace_line(folder + "/ToyConfiguration.txt", 17, "CosmicRate = %f\n" % (rate[i]))
		replace_line(folder + "/ToyConfiguration.txt", 18, "Efficiency = %f\n" % (efficiency[i]))

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
trial = 999
Nsum = 5

def task(i):
	values = ROOT.ScanAnalysis(folder + str(i) + "/", folder +  str(i) + "/ToyConfiguration.txt",0,trial,3,0.3, Nsum,rate[i])
	npResult = np.asarray(values)
	return npResult

if __name__ == "__main__":
	# create the process pool
	index = range(0,len(rate))
	with Pool(processes = 10, maxtasksperchild = 1) as pool:
		for results in pool.map(task, index):
			#print(results)
			bias_thr = np.append(bias_thr, results[0])
			variance_thr = np.append(variance_thr, results[1])
			bias_fw = np.append(bias_fw, results[2])
			variance_fw = np.append(variance_fw, results[3])
			bias_rev = np.append(bias_rev, results[4])
			variance_rev = np.append(variance_rev, results[5])
			bias_cf = np.append(bias_cf, results[6])
			variance_cf = np.append(variance_cf, results[7])
			bias_sum = np.append(bias_sum, results[8])
			variance_sum = np.append(variance_sum, results[9])


# define a cost function, function of bias and variance
cost_cf = (np.abs(bias_cf) + variance_cf)**2/variance_cf	
cost_fw = (np.abs(bias_fw) + variance_fw)**2/variance_fw
cost_rev = (np.abs(bias_rev) + variance_rev)**2/variance_rev
cost_sum = (np.abs(bias_sum) + variance_sum)**2/variance_sum
cost_thr = (np.abs(bias_thr) + variance_thr)**2/variance_thr

plt.figure(1, figsize = (12,12))
plt.title("MVA scan " + r"$N_{trial}$ = " + str(trial + 1) , fontsize = 20)
plt.xlabel("rate [events/second]", fontsize = 15)
plt.ylabel("cost [kHz]", fontsize = 15)
plt.xscale("log")
plt.grid()
plt.errorbar(rate,cost_cf, linestyle = '',marker = "v", color = "red", label = "constant fraction 30%")
plt.errorbar(rate,cost_fw, linestyle = '',marker = ".", color = "blue", label = "foward")
plt.errorbar(rate,cost_rev, linestyle = '',marker = "D", color = "green", label = "reversed")
plt.errorbar(rate,cost_sum, linestyle = '',marker = "s", color = "black", label = "Neighbors sum (N = %d)" % Nsum)
plt.errorbar(rate,cost_thr, linestyle = '',marker = "D", color = "orange", label = "over Threshold (> 3)")
plt.legend(fontsize = 10)

plt.figure(2, figsize = (12,12))
plt.title("Bias " + r"$N_{trial}$ = " + str(trial + 1) , fontsize = 20)
plt.xlabel("rate [events/second]", fontsize = 15)
plt.ylabel("Bias", fontsize = 15)
plt.xscale("log")
plt.grid()
plt.errorbar(rate,np.abs(bias_cf), linestyle = '',marker = "v", color = "red", label = "constant fraction 30%")
plt.errorbar(rate,np.abs(bias_fw), linestyle = '',marker = ".", color = "blue", label = "foward")
plt.errorbar(rate,np.abs(bias_rev), linestyle = '',marker = "D", color = "green", label = "reversed")
plt.errorbar(rate,np.abs(bias_sum), linestyle = '',marker = "s", color = "black", label = "Neighbors sum (N = %d)" % Nsum)
plt.errorbar(rate,np.abs(bias_thr), linestyle = '',marker = "D", color = "orange", label = "over Threshold (> 3)")
plt.legend(fontsize = 10)

plt.figure(3, figsize = (12,12))
plt.title(r"$\sigma$ vs Bias " + r"$N_{trial}$ = " + str(trial + 1) , fontsize = 20)
plt.xlabel(r"$\sigma$ [kHz]", fontsize = 15)
plt.ylabel("Bias [kHz]", fontsize = 15)
plt.grid()
plt.errorbar(variance_cf,bias_cf, linestyle = '',marker = "v", color = "red", label = "constant fraction 30%")
plt.errorbar(variance_fw,bias_fw, linestyle = '',marker = ".", color = "blue", label = "foward")
plt.errorbar(variance_rev,bias_rev, linestyle = '',marker = "D", color = "green", label = "reversed")
plt.errorbar(variance_sum,bias_sum, linestyle = '',marker = "s", color = "black", label = "Neighbors sum (N = %d)" % Nsum)
plt.errorbar(variance_thr,bias_thr, linestyle = '',marker = "D", color = "orange", label = "over Threshold (> 3)")
plt.legend(fontsize = 10)

plt.figure(4, figsize = (12,12))
plt.title("MVA scan " + r"$N_{trial}$ = " + str(trial + 1) , fontsize = 20)
plt.xlabel("rate [events/second]", fontsize = 15)
plt.ylabel(r"$\sigma$ [kHz]", fontsize = 15)
plt.xscale("log")
plt.grid()
plt.errorbar(rate,variance_cf, linestyle = '',marker = "v", color = "red", label = "constant fraction 30%")
plt.errorbar(rate,variance_fw, linestyle = '',marker = ".", color = "blue", label = "foward")
plt.errorbar(rate, variance_rev, linestyle = '',marker = "D", color = "green", label = "reversed")
plt.errorbar(rate,variance_sum, linestyle = '',marker = "s", color = "black", label = "Neighbors sum (N = %d)" % Nsum)
plt.errorbar(rate,variance_thr, linestyle = '',marker = "D", color = "orange", label = "over Threshold (> 3)")
plt.legend(fontsize = 10)
plt.show()
