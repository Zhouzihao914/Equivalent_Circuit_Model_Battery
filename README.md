# Simplified_ECM_Matlab

Recreate the **Thevenin model** for a battery. Neglect the temperature influence. \
Use the data set from Prof. Gregory L Plett's book <[Battery Modelling](http://mocha-java.uccs.edu/BMS1/index.html)>.\
![Image of ECM](https://gitlab.com/zihaos-play-yard/simplified_ecm_matlab/-/blob/main/ECM.PNG)

## Getting started

Begin with the onesampleOCV.m file to get required static parameters (R0, OCV-SOC, eta, Q_tol).\
Then run the onesampleDynamic.m file to get required dynamic parameters (OCV, eta, R0, RC).\
Finally test the built ECM model with testECM.m file.\
ECMcell.m -- Simulates an ECM model for input current ik at temperature T.\
AgingFuncdisQ.m -- define the discharge capacity into a function of time. \
simCCCV.m /simCPCV.m -- simulate Constant Current Constant Voltage(CC-CV) and Constant Power Constant Voltage(CP-CV) charging profile.\
generateSynthetic.m -- generating synthetic data which simulate the real usage of battery in Electrical vehicles.\

## Project status
Most parts are done. May restruct some codes later for object oriented programming.\

## Authors and acknowledgment
Author: Zihao Zhou




