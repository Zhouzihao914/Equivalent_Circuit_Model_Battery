# Simplified_ECM_Matlab

Recreate the **Thevenin model** for a battery. Neglect the temperature influence. \
Use the data set from Prof. Gregory L Plett's book <[Battery Modelling](http://mocha-java.uccs.edu/BMS1/index.html)>.\
![Image of ECM](https://gitlab.com/zihaos-play-yard/simplified_ecm_matlab/-/blob/main/ECM.PNG)

## Getting started

Begin with the onesampleOCV.m file to get required static parameters (R0, OCV-SOC, eta, Q_tol).
Then run the onesampleDynamic.m file to get required dynamic parameters (OCV, eta, R0, RC).
Finally test the built ECM model with testECM.m file.

There will be a AgingFuncdisQ.m file to define the discharge capacity into a function of time. 
There will be simulation file for CC-CV and CP-CV profiles.
There will be generateSynthetic.m for generating synthetic data which simulate the real usage of battery in Electrical vehicles.

## Project status
To be continued...

## Authors and acknowledgment
Author: Zihao Zhou




