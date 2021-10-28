# Simplified_ECM_Matlab

Recreate the **Thevenin model** for a battery. Neglect the temperature influence. 
Use the data set from Prof. Gregory L Plett's book <[Battery Modelling](http://mocha-java.uccs.edu/BMS1/index.html)>.
![Image of ECM](https://gitlab.com/zihaos-play-yard/simplified_ecm_matlab/-/blob/main/ECM.PNG)

## Getting started

Begin with the onesampleOCV.m file to generate required static parameters (Resistance, OCV-SOC, eta, Q_tol).
Then run the onesampleDynamic.m file to train required dynamic parameters (I_R1, SOC).
Finally test the built ECM model with testModel.m file.

There will be a tol_Q_func.m file to define the total capacity into a function of time. 

## Project status
Daily Update. To be continued...

## Authors and acknowledgment
Author: Zihao Zhou




