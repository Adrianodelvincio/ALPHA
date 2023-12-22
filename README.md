# ALPHA
### Software Code of Monte Carlo Simulation for ALPHA-2 Hyperfine measurement
#### PhD student Adriano Del Vincio (University of Brescia)


This repository contains software developed for ALPHA-2 and ALPHA-g experiments. We have developed a model for distinguishig annihilations of Trap Walls from residual gas annihilation.

Language code: **ROOT**, **Pyroot** and **Python**

Structure of the repository
---
**ALPHA-2** folder contains both the data and the analysis code required to build and evaluate a new model for the radius distribution of annihilation events in ALPHA experiment. This new model is designed to account for the distinct contributions of the annihilation on trap walls, annihilation due to residual gas, and background from muons. 

1. The data that are used in the analysis are in **.root** format and are saved in the folder **DataSetROOT**. The datasets analyzed include an almost pure sample of annihilation on walls (*..._mixing.vertex.root*), an almost pure sample of annihilation due to residual gas (*..._uwlosses_160.vertex.root*) and datasets cointaing only cosmic background events (*.._cosmic.vertex.root*) 



---

The software is dived in different scripts, the important ones for building up a model to distinguish the different mechanism of annihilation are

1. *TemplateMLE.cpp*
2. *AnalyticMLE.cpp*
2. *MonteCarloForFit.cpp*
3. *Conversion.cpp*
4. *fitAllFrequency.cpp*

1) In *TemplateMLE.cpp* an MLE fit is performed to the *cosmic*, *uwlosses* and *mixing* datasets. The data are analyzed with RDataFrame. The radius of each annihilation events is computed from the X and Y vertex coordinates. Analytic models are used to fit the radius distributions obatined for the anti-hydrogen annihilations on walls/residual gas and for the cosmic events.

2) In *AnalyticMLE.cpp* the analytic models defined in *TemplateMLE.cpp* are combined together. This final model can be used to analyze runs which have different percentages of Walls/residual gas annihilations.

The program can be run passing as argument the name of the run to be analyzed.

root [0] .x AnalyticMLE.cpp("filename")

3) In *MonteCarloForFit.cpp* we have developed a simulation to study the fit procedure applied in *AnalyticMLE.cpp*. The program generates fake data that can be used to validate the fit procedure. The proportion of simulated data resulting from residual gas/wall annihilation is regulated by a parameter of the program. The analytic model is used to fit the simulated data and retrieve the number of annihilation events belonging to the different process.

The parameters of the simulation are loaded from a separated *txt* files, named *configToyModel.txt*, which contains:

1. N: number of annihilation per run
2. Nloop: number of run which are generated and fitted
3. a: percentage of annihilation on walls
4. Ncomsic: number of cosmic events
5. mu: parameter of the analytic model for the annihilation on walls
6. sigWall: parameter of the analytic model for the annihilation on walls
7. sigRay: parameter of the analytic model for the annihilation from res. gas



---
**LineShape** simulation

the folder **LineShape** contains the code for a Monte Carlo toy of the hyperfine measurements. The aim of the simulation is to test different onset finding algorithm, studying the bias and statistical uncertainty. The main structure of the simulation is presented here:
![alt text](https://github.com/Adrianodelvincio/ALPHA/blob/main/ALPHA-2/LineShape/Beamer/SimulationScheme.pdf)
