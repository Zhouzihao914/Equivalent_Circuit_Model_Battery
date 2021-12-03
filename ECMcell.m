% Simulates an ECM model for input current ik at temperature T
%
% Inputs:  ik - current, where (+) is discharge
%          T  - temperature (degC)
%      deltaT - sampling interval in data (s)
%       model - standard model structure
%          z0 - initial cell state of charge
%         iR0 - initial resistor currents as column vector
%          Q0 -
%           N - cyclenum 
% Outputs: vk - cell voltage for all timesteps
%         irk - resistor currents (in R-C branches) for all timesteps
%          hk - hysteresis states for all timesteps
%          zk - state of charge for all timesteps
%         sik - instantaneous hysteresis for all timesteps
%         OCV - open-circuit voltage for all timesteps

function [vk,irk,zk,OCV,Qend] = ECMcell(ik,T,deltaT,model,z0,iR0,Q0,func_ind,a,b,N)
  format long
  % Force data to be column vector(s)
  ik = ik(:); iR0 = iR0(:);
  % Get model parameters from model structure
  RCfact = exp(-deltaT./abs(getParamESC('RCParam',T,model)))';
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
  
  % Assume the ECM model have built in toltal capacity time func.
  Ltime = length(ik);
  delta_Q = AgingFuncDisQ(func_ind, a, b, N, Ltime);
  if isscalar(delta_Q)
    delta_Q = ones([Ltime,1])*delta_Q;
  end
  Q_time = 3600*(Q0 - delta_Q);
  zk = z0-cumsum([0;etaik(1:end-1)])*deltaT./Q_time;
  Qend = Q_time/3600;

  if any(zk>1.1)
    warning('Current may have wrong sign as SOC > 110%');
  end
  
  if any(zk<-.1)
    warning('Current may have wrong sign as SOC < -10%');
  end
    
  % Compute output equation
  OCV = getOCVfromSOCTemp(zk,T,model);
  vk = OCV - irk*RParam' - ik.*R0Param;

