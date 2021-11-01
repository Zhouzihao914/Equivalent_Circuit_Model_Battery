% Simulates an ECM model for input current ik at temperature T
%
% Inputs:  ik - current, where (+) is discharge
%          T  - temperature (degC)
%      deltaT - sampling interval in data (s)
%       model - standard model structure
%          z0 - initial cell state of charge
%         iR0 - initial resistor currents as column vector
%          h0 - initial hysteresis state
% Outputs: vk - cell voltage for all timesteps
%         irk - resistor currents (in R-C branches) for all timesteps
%          hk - hysteresis states for all timesteps
%          zk - state of charge for all timesteps
%         sik - instantaneous hysteresis for all timesteps
%         OCV - open-circuit voltage for all timesteps


function [vk,irk,zk,OCV] = ECMcell(ik,T,deltaT,model,z0,iR0)
  % Force data to be column vector(s)
  ik = ik(:); iR0 = iR0(:);
  % Get model parameters from model structure
  RCfact = exp(-deltaT./abs(getParamESC('RCParam',T,model)))';
  G = getParamESC('GParam',T,model);
  Q = getParamESC('QParam',T,model);
  RParam = getParamESC('RParam',T,model);
  R0Param = getParamESC('R0Param',T,model);
  etaParam = getParamESC('etaParam',T,model);
  
  etaik = ik; 
  etaik(ik<0) = etaParam*ik(ik<0); % compensate for coulombic efficiency
  
  % Simulate the dynamic states of the model
  if exist('ss','file'), % use control-system-toolbox method, if available
    sysd= ss(diag(RCfact),1-RCfact,eye(length(RCfact)),0,-1);
    irk = lsim(sysd,etaik,[],iR0);
  else
    irk=zeros([length(ik) length(iR0)]); irk(1,:) = iR0;
    for k = 2:length(ik),
      irk(k,:) = RCfact'.*irk(k-1,:) + (1-RCfact')*etaik(k-1);
    end
  end
  zk = z0-cumsum([0;etaik(1:end-1)])*deltaT/(Q*3600); 
  if any(zk>1.1),
    warning('Current may have wrong sign as SOC > 110%');
  end
    
  % Compute output equation
  OCV = getOCVfromSOCTemp(zk,T,model);
  vk = OCV - irk*RParam' - ik.*R0Param;