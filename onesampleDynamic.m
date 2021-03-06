% script onesampleDynamic.m
% load the data from a cell <A123> tested under eight different 
% temperatures. This script calls generateDynamic.m to calculate dynamic 
% parameters (R0/RC), then saves the model to a model file. Please run
% onesampleOCV.m before this file.

% This script is a self-reproduced simplified version from the one provided
% by Prof. Plett, Gregory L. <http://mocha-java.uccs.edu/BMS1/index.html>

close all
clear

cellID = 'A123';
num_rc = 1;
temps = [-25  -15   -5    5   15   25   35   45];
I_mags =  [10   10   30   45   45   50   50   50];

modelFile = sprintf('%smodel-ocv.mat',cellID);
if ~exist(modelFile, 'file'),
    error(['The OCV-SOC model file "%s" does not exist in current path.\n'...
        'Please change into the right folder.'],modelFile);
end
load(modelFile);
% I_mags represent the maximum C-rate which should be positive
data = zeros([0 length(I_mags > 0)]); 
data_ind = 0;
  for temps_ind = 1:length(temps),
    curr_mag = I_mags(temps_ind);     
    if curr_mag < 0,                     
      continue 
    else                                
      data_ind = data_ind + 1;
    end
    % if temperature is negative, then load this
    if temps(temps_ind) < 0, 
      DYNPrefix = sprintf('%s_DYN/%s_DYN_%02d_N%02d',cellID,cellID,...
                          curr_mag,abs(temps(temps_ind)));
    else
    % if temperature is positive, then load this
      DYNPrefix = sprintf('%s_DYN/%s_DYN_%02d_P%02d',cellID,cellID,...
                          curr_mag,temps(temps_ind));
    end
    inFile = sprintf('%s.mat',DYNPrefix);
    if ~exist(inFile,'file'),
      error(['File "%s" not found.\n' ...
        'Please change folders so that "%s" points to a valid data '...
        'file and re-run runProcessDynamic.'],inFile,inFile); %#ok<SPERR>
    end
    fprintf('Loading %s\n',inFile); load(inFile);        
    data(data_ind).temp    = temps(temps_ind); % store temperature
    data(data_ind).script1 = DYNData.script1; % store data from each of the
    data(data_ind).script2 = DYNData.script2; % three scripts
    data(data_ind).script3 = DYNData.script3;
  end

  model = generateDynamic(data,model,num_rc);
  save(modelFile,'model');

  % Plot model-match voltage results at 25 degC, plus RMS voltage-estimation 
  % error between 5% and 95% cell state of charge
  ind_25 = find(temps == 25);
  [vk,rck,zk,OCV] = ECMcell(data(ind_25).script1.current,temps(ind_25),1,...
                            model,1,zeros(num_rc,1),model.QParam(ind_25),0,0,0);

  figure(4);
  tk = (1:length(data(ind_25).script1.current))-1;
  plot(tk/60,data(ind_25).script1.voltage,tk/60,vk);
  xlabel('Time(min)')
  ylabel('Voltage(volts)')
  title('Voltage and estimates at T=25');
  legend('real voltage','estimated voltage');

  tol_verr = [];
  for i = 1:length(temps),
    [vk,rck,zk,OCV] = ECMcell(data(i).script1.current,temps(i),1,...
                            model,1,zeros(num_rc,1),model.QParam(ind_25),0,0,0);
    curr_verr = data(i).script1.voltage - vk';
    v1 = getOCVfromSOCTemp(0.95,temps(i),model);
    v2 = getOCVfromSOCTemp(0.05,temps(i),model);
    N1 = find(data(i).script1.voltage<v1,1,'first'); 
    N2 = find(data(i).script1.voltage<v2,1,'first');
    if isempty(N1), N1=1; end; if isempty(N2), N2=length(curr_verr); end
    curr_mse=sqrt(mean(curr_verr(N1:N2).^2));
    tol_verr = [tol_verr, curr_mse];
  end
  tol_verr = tol_verr*1000;
  disp('RMS error of ECMcell at:');
  disp(temps);
  disp(tol_verr);

  figure(5)
  plot(temps,tol_verr,'ro');
  xlabel('Temperature');
  ylabel('Mean Square Error');

  figure(6)
  plot(tk/60,data(ind_25).script1.voltage);
  xlabel('Time(min)')
  ylabel('Voltage(volts)')
  title('Terminal Voltage at T=25');


