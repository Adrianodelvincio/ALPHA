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

1) In *TemplateMLE.cpp* an MLE fit is performed to the *cosmic*, *uwlosses* and *mixing* datasets. Some analytic functions are used to fit the radius distribution of the anti-hydrogen annihilations on walls/residual gas and cosmic events. These models are used in combination to fit the data of the annihilation contained in *DataSetROOt* folder. In the end we extract from the data a template model of the radius distribution. The templates are combined together to build a rough model which can be used to fit the data.

2) In *AnalyticMLE.cpp* the analytic models used to fit the radius distributions in *TemplateMLE.cpp* are combine together to build a global analytic model to fit the data.

3) To study the feasibility of this method, we have developed a simple Monte Carlo simulation where we generates some data which follows the radius distributions analysed in *TemplateMLE.cpp* and we apply our fit procedure. We study the uncertainty and bias associated of the fit model. 


---
**LineShape** simulation

the folder **LineShape** contains the code for a Monte Carlo toy of the hyperfine measurements. The aim of the simulation is to test different onset finding algorithm, studying the bias and statistical uncertainty. The main structure of the simulation is presented here:
![alt text](https://github.com/Adrianodelvincio/ALPHA/blob/main/ALPHA-2/LineShape/Beamer/SimulationScheme.pdf)
