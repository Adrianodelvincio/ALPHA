DATA TYPE IN DATASET

type 0 -> annihilation on wall
type 1 -> annihilation on residual gas
type 2 -> cosmic background


To run the analysis, call ROOT and run the script Analysis.cpp

For instance, if you want to run the analysis on a set of data, as the first 100 dataset contained in the directory "folder", run: 

root [1] .x  Analysis.cpp("folder/",  0,  100,  10,  0.3, 4, 5, 5, 5,"Plot/")

1) The folder where the data are saved
2) starting data to be analyzed
3) last data be analyzed
4) Nfilter parameter
5) Fraction parameter for constant fraction algorithm
6) Threshold for signle threshold algorithm
7) First threshold for forward or reversed algorithm
8) Second threshold for forward or reversed algorithm
8) Path to the folder where to save the plots

!ACHTUNG!
Check if the datafolder contains the toyconfiguration files, needed to load the parameters for the analysis!

### LIST OF AVAILABLE DATASET SAVED

- scanMvaData

this folder contains all the datasets generated for each MVA point. The hyperfine splitting is uniformly generated around the cpt value with interval (-30, + 10) kHz.

- scanMvaData2

this folder contains all the datasets generated for each MVA point. The hyperfine splitting is uniformly generated around the cpt value with interval (-15, + 15) kHz.

- point 1 

Data corresponding to MVA working point rate = 0.046, efficiency = 0.775

- point 2

Data corresponding to MVA working point rate = 0.01, efficiency = 0.638

- point 3

Data corresponding to MVA working point rate = 0.23, efficiency = 0.86

- Test_1

Test 1 of the cost function.

Interval of generation:  (-20,20)

- Test_2

Test 2 of the cost function

interval of generation [-40,0] (pay attention, binbeforeonset = 4)

- Test_3

Test 3 of the cost function

interval of generation [-40, 40] (pay attention, binbeforeonset = 4)

- Test_4 

Test 4 of the cost function

interval of generation [-40, 20] (pay attention, binbeforeonset = 4)

------------------------------------------------------------------------
------------------------------------------------------------------------
FROM HERE, ALL THE DATASET HAVE BDRIFT

- Passcut

Passcut dataset, The hyperfine splitting is uniformly generated around the cpt value with interval (-30, + 10) kHz.


- scanMvaData3

this folder contains all the datasets generated for each MVA point. The hyperfine splitting is uniformly generated around the cpt value with interval (-30, + 10) kHz. The Bdrift effects are simulated, also the start of the sweep is considered with operator error.

-ScanMvaData4

Same ad scanMvaData3 but ALL the Bdrift effects are simulated.

- longWindow

Bdrift put equal to zero, sweep step equal to 100 (long windows)
