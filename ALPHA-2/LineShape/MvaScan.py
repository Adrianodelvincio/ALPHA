import ROOT
import numpy as np
import array
import os

rate, efficiency = np.loadtxt("MVA-points.txt", unpack = True)

# DATA GENERATION
simulationCPP = "LoopLineShape"
mvaScan = "\"true\""
Nfiles = "500"
ConfFile = "\"ToyConfiguration.txt\""

ROOT.gInterpreter.ProcessLine(".L LoopLineShape.cpp")
"""
for i, item in enumerate(efficiency):
	CosmicRate = str(rate[i])
	Efficiency = str(item)
	folder = "mva_" + str(i)
	# Create folder
	
	if not os.path.exists("/home/adriano/Documents/ALPHA/ALPHA-2/LineShape/" + folder):	# Check if the folder exist
		os.makedirs("/home/adriano/Documents/ALPHA/ALPHA-2/LineShape/" + folder)			# Create the folder

	code = simulationCPP + "(" + Nfiles + ",\"" + folder + "/\"," + mvaScan + "," + CosmicRate + "," + Efficiency + "," + ConfFile + ")"
	print("running: ", code)
	# Execute simulation code
	ROOT.gInterpreter.ProcessLine(code)
	# copy the configuration file in the folder
	os.popen("cp ToyConfiguration.txt " + folder + "/"  + ConfFile)
"""
for i in range(18,len(efficiency)):
	CosmicRate = str(rate[i])
	Efficiency = str(efficiency[i])
	folder = "mva_" + str(i)
	# Create folder
	
	if not os.path.exists("/home/adriano/Documents/ALPHA/ALPHA-2/LineShape/" + folder):	# Check if the folder exist
		os.makedirs("/home/adriano/Documents/ALPHA/ALPHA-2/LineShape/" + folder)			# Create the folder

	code = simulationCPP + "(" + Nfiles + ",\"" + folder + "/\"," + mvaScan + "," + CosmicRate + "," + Efficiency + "," + ConfFile + ")"
	print("running: ", code)
	# Execute simulation code
	ROOT.gInterpreter.ProcessLine(code)
	# copy the configuration file in the folder
	os.popen("cp ToyConfiguration.txt " + folder + "/"  + ConfFile)	

# DATA ANALYSIS
# Execute analysis code
"""
ROOT.gInterpreter.ProcessLine(".L ScanAnalysis.cpp")
for i, item in enumerate(rate):
	analysisCPP = "ScanAnalysis"
	folder = "\"mva_" + str(i) + "/\""
	ConfFile = "\"mva_" + str(i) + "/ToyConfiguration.txt\""
	listFile = "0,10"
	algorithmParameters = "1,0.1," + str(rate[0])
	#result = ROOT.std.vector("double")(10)
	result = ROOT.gInterpreter.ProcessLine(analysisCPP + "(" + folder + "," + ConfFile + "," + listFile + "," + algorithmParameters + ")")
	print("point %d" % i, " ", type(result), " ", result)
"""
