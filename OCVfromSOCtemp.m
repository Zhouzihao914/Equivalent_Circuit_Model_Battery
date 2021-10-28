% function ocv = OCVfromSOCtemp(soc,temp,model)
%
% Computes the fully rested open-circuit voltage for a particular state of 
% charge and temperature point
%
% Inputs: soc = scalar or matrix of cell state of charge points
%        temp = scalar or matrix of temperatures at which to calc. OCV
%       model = data structure produced by processDynamic
% Output: ocv = scalar or matrix of open circuit voltages -- one for every
%               soc and temperature input point

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

function ocv=OCVfromSOCtemp(soc,temp,model)
  soccol = soc(:); % force soc to be column vector
  SOC = model.SOC(:); % force to be column vector
  OCV0 = model.OCV0(:); % force to be column vector
  OCVrel = model.OCVrel(:); % force to be column vector
  if isscalar(temp), % replicate a scalar temperature for all socs
    tempcol = temp*ones(size(soccol)); 
  else
    tempcol = temp(:); % force matrix temperature to be column vector
    if ~isequal(size(tempcol),size(soccol)),
      error(['Function inputs "soc" and "temp" must either have same '...
        'number of elements, or "temp" must be a scalar']);
    end
  end
  diffSOC=SOC(2)-SOC(1); % spacing between SOC points - assume uniform
  ocv=zeros(size(soccol)); % initialize output to zero
  I1=find(soccol <= SOC(1)); % indices of socs below model-stored data
  I2=find(soccol >= SOC(end)); % and of socs above model-stored data
  I3=find(soccol > SOC(1) & soccol < SOC(end)); % the rest of them
  I6=isnan(soccol); % if input is "not a number" for any locations

  % for voltages less than lowest stored soc datapoint, extrapolate off 
  % low end of table 
  if ~isempty(I1),
    dv = (OCV0(2)+tempcol.*OCVrel(2)) - (OCV0(1)+tempcol.*OCVrel(1));
    ocv(I1)= (soccol(I1)-SOC(1)).*dv(I1)/diffSOC + ...
             OCV0(1)+tempcol(I1).*OCVrel(1);
  end

  % for voltages greater than highest stored soc datapoint, extrapolate off
  % high end of table
  if ~isempty(I2),
    dv = (OCV0(end)+tempcol.*OCVrel(end)) - ...
         (OCV0(end-1)+tempcol.*OCVrel(end-1));
    ocv(I2) = (soccol(I2)-SOC(end)).*dv(I2)/diffSOC + ...
              OCV0(end)+tempcol(I2).*OCVrel(end);
  end

  % for normal soc range, manually interpolate (10x faster than "interp1")
  I4=(soccol(I3)-SOC(1))/diffSOC; % using linear interpolation
  I5=floor(I4); I45 = I4-I5; omI45 = 1-I45;
  ocv(I3)=OCV0(I5+1).*omI45 + OCV0(I5+2).*I45;
  ocv(I3)=ocv(I3) + tempcol(I3).*(OCVrel(I5+1).*omI45 + OCVrel(I5+2).*I45);
  ocv(I6)=0; % replace NaN SOCs with zero voltage
  ocv = reshape(ocv,size(soc)); % output is same shape as input