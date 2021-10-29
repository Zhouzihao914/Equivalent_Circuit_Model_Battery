
% Function getOCVfromSOCTemp
% Basically, it is a lookup table for OCV from customize socs and temps.
% So, it requires the built reference OCV-SOC relationship.

% Inputs:
%   socs = the soc array 
%   temps = the temperature array
%   model = the trained model, specifically, the OCV-SOC
%   relationship trained from generateOCVSOC
% Output:
%   OCV values corresponding to socs and temps

function ocv = getOCVfromSOCTemp(socs, temps, model)
    soc_col = socs(:);
    SOC = model.SOC(:);
    OCV0 = model.OCV0(:); 
    OCVrel = model.OCVrel(:);
    
    % temps may be a single scalar which represent constant temperature
    % also it may be a array which represent variable temperature
    if isscalar(temps),
        temps_col = temps*ones(size(soc_col));
    else
        temps_col = temps(:);
        if ~isequal(size(temps_col),(soc_col)),
            error(['Function inputs "socs" and "temps" must either have same '...
            'number of elements, or "temp" must be a scalar'])
        end
    end
    
    % Assume the soc interval is uniform
    % Divide the array into different portions according to values
    soc_interval = SOC(2) - SOC(1);
    ocv = zeros(size(soc_col));
    ind_tail = find(soc_col <= SOC(1));
    ind_hat = find(soc_col >= SOC(end));
    ind_mid = find(soc_col > SOC(1) & soc_col < SOC(end));
    ind_nan = isnan(soc_col);
    
    % For soc levels below the recorded minimum SOC, extrapolate off low end of
    % the recorded SOC.
    if ~isempty(ind_tail),
        slope_tail = ((OCV0(2) + temps_col.*OCVrel(2)) - ...
                     (OCV0(2) - temps_col.*OCVrel(2)))/soc_interval;
        ocv(ind_tail) = (soc_col(ind_tail) - SOC(1)).* slope_tail(ind_tail) + ...
                        OCV0(1) + temps_col(ind_tail).*OCVrel(1);
    end
    
    % For soc levels above the recorded maximum SOC, extrapolate off high end of
    % the recorded SOC.
    if ~isempty(ind_hat),
        slope_hat = ((OCV0(end) + temps_col.*OCVrel(end)) - ...
                    (OCV0(end-1) + temps_col.*OCVrel(end-1)))/soc_interval;
        ocv(ind_hat) = (soc_col(ind_hat) - SOC(end)).* slope_hat(ind_hat) + ...
                       OCV0(end) + temps_col(ind_hat).*OCVrel(end);
    end
    
    % For normal soc range, manually interpolate (10x faster than "interp1")
    % The intuition is that the actual soc value may between recorded SOC
    % points, such as 0.5%(record SOC)<0.7%(soc value)<1.0%(record SOC). so ocv
    % value should be a linear combination from the floor and ceil record OCV.
    % OCV(soc value) = (1-c)*OCV(floor record SOC) + c*OCV(ceil record SOC
    % c = distance(soc value - floor record SOC) - [0,1]
    
    soc_actual_ratio = (soc_col(ind_mid) - SOC(1))/soc_interval;
    soc_record_floor = floor(soc_actual_ratio);
    distance_floor = soc_actual_ratio - soc_record_floor;
    ocv(ind_mid) = OCV0(soc_record_floor + 1).*(1-distance_floor) + ...
                   OCV0(soc_record_floor + 2).*distance_floor;
    ocv(ind_mid) = ocv(ind_mid) + temps_col(ind_mid).*( ...
                   OCVrel(soc_record_floor + 1).*(1-distance_floor) + ...
                   OCVrel(soc_record_floor + 2).*distance_floor);
    % For NaN soc values, replace it with zero
    ocv(ind_nan) = 0;
    ocv = reshape(ocv, size(soc_col));












