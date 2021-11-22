% Simulates a Constant Current Constant Voltage charging profile for an
% established ECM cell. Running generateOCVSOC.m and generateDynamic.m
% before running this script

clear; close all;

modelFile = 'A123model-ocv.mat';
load(modelFile);

%Get or set parameters for the ECM model
maxtime = 1501; T = 25;
q = getParamESC('QParam',T,model);
rc = exp(-1./abs(getParamESC('RCParam',T,model)));
r  = (getParamESC('RParam',T,model));
r0 = getParamESC('R0Param',T,model);
maxV = 3.5;

storez = zeros([maxtime 1]);  % storage for SOC
storev = zeros([maxtime 1]);  % storage for voltage
storei = zeros([maxtime 1]);  % storage for current
storep = zeros([maxtime 1]);  % storage for power
z  = 0.1; irc = zeros([size(model.RCParam,2),1]); % initialize to 50% SOC, resting
CC = 10;

% Simulate CC/CV
for k = 1:maxtime,
  v = getOCVfromSOCTemp(z,T,model) - r*irc; % fixed voltage

  ik = (v - maxV)/r0; % compute test ik to achieve maxV
  ik = max(-CC,ik);   % but limit current to no more than CC in mag.

  z = z - (1/3600)*ik/q;  % Update cell SOC
  irc = irc.*rc' + (1-rc)'.*ik; % Update resistor currents

  storez(k) = z; % Store SOC for later plotting
  storev(k) = v - ik*r0;
  storei(k) = ik; % store current for later plotting
  storep(k) = ik*storev(k);
end % for k

time = 0:maxtime -1;
subplot(2,2,1); plot(time,100*storez); 
title('State of charge versus time');
xlabel('Time (s)'); ylabel('SOC (%)'); grid on

subplot(2,2,2); plot(time,storev); 
title('Terminal voltage versus time');
xlabel('Time (s)'); ylabel('Voltage (V)');
grid on

subplot(2,2,3); plot(time,storei); 
title('Cell current versus time');
xlabel('Time (s)'); ylabel('Current (A)');
grid on

subplot(2,2,4); plot(time,storep);
title('Cell power versus time');
xlabel('Time (s)'); ylabel('Power (W)');
grid on
