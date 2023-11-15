#ALPHA
### Project code for analysis of ALPHA-g experiment
#### PhD student Adriano Del Vincio


This repository contains the software for the analysis of the ALPHA-g experiment.

Language code: **ROOT**, **Pyroot** and **Python**

Structure of the repository
---
**ALPHA-2** folder contains the data and the analysis code for the ALPHA-2 spectroscopy experiment. 

1. the data are a series of.cvs files. The files that contains X,Y,Z position of the annihilation vertices are the ones that contains the word **vertex**. For instance *r69177_cosmics.vertex.csv* contains the data about the vertices of cosmic background
2. The data files are saved in **Control** folder and **Dataset** folder, that are not saved in this repository. In the control folder we have two subfolder: **cosmic** and **intrap**. In **cosmic**, as the word says, the data are used to study the distribution of the cosmic background of the experiment (mainly muons). **intrap** contains the data about the annihilation due to the residual gas (indicated as **..._uwlosses_160.vertex.csv**) and almost pure data of annihilation on the trap walls, indicated as **..._mixing.vertex.csv**
---

The software is dived in different scripts, that are

1. *AnalysisMLE.cpp*
2. *MonteCarloForFit.cpp*
3. *Conversion.cpp*
4. *fitAllFrequency.cpp*
5. *GaussianToy.cpp*
