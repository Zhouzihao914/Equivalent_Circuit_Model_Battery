% Function generateOCVSOC
% Inputs:
%   data = cell-test data passed in from onesampleOCV
%   cellID = cell identifier (string)
%   minV = minimum cell voltage to use in OCV relationship
%   maxV = maximum cell voltage to use in OCV relationship
% Output:
%   model = data structure with information for recreating OCV

function model = generateOCVSOC(data, cellID, minV, maxV)
    filetemps = [data.temp]; filetemps = filetemps(:);
    numtemps = length(filetemps); 

    % Use data under 25 degC to calculate coulombic efficiency and capacity
    ind25 = find(filetemps == 25); 
    if isempty(ind25),
        error('data under 25 degC is required!')
    end
    indnot25 = find(filetemps~=25);
    SOC = 0:0.005:1;
    file_data = zeros([0 length(data)]); % record test data and processed data
    eta = zeros([0 length(data)]);
    tol_Q = zeros([0 length([0 length(data)])]);
    tol_disAh = data(ind25).script1.disAh(end) + data(ind25).script2.disAh(end) + ...
                data(ind25).script3.disAh(end) + data(ind25).script4.disAh(end);
    tol_ChgAh = data(ind25).script1.chgAh(end) + data(ind25).script2.chgAh(end) + ...
                data(ind25).script3.chgAh(end) + data(ind25).script4.chgAh(end);

    % 25 degC coulombic efficiency
    eta_25 = tol_disAh/tol_ChgAh; eta(ind25) = eta_25;
    % adjust charge Ah in all scripts per eta25
    data(ind25).script1.chgAh = data(ind25).script1.chgAh*eta_25;
    data(ind25).script2.chgAh = data(ind25).script2.chgAh*eta_25; 
    data(ind25).script3.chgAh = data(ind25).script3.chgAh*eta_25; 
    data(ind25).script4.chgAh = data(ind25).script4.chgAh*eta_25;
    
    % 25 degC total capacity
    Q_25 = data(ind25).script1.disAh(end) + data(ind25).script2.disAh(end) - ...
            data(ind25).script1.chgAh(end) - data(ind25).script2.chgAh(end);
    tol_Q(ind25) = Q_25;
    
    % Resistants from both beginning(R1) and end(R2) of slow discharge step
    ind_dis  = find(data(ind25).script1.step == 2);
    IR1_dis = data(ind25).script1.voltage(ind_dis(1)-1) - ... 
              data(ind25).script1.voltage(ind_dis(1)); 
    IR2_dis = data(ind25).script1.voltage(ind_dis(end)+1) - ... 
              data(ind25).script1.voltage(ind_dis(end)); 

    % Resistants from both beginning(R1) and end(R2) of slow charge step
    ind_chg = find(data(ind25).script3.step == 2);
    IR1_chg = data(ind25).script3.voltage(ind_chg(1)) - ... 
              data(ind25).script3.voltage(ind_chg(1)-1); 
    IR2_chg = data(ind25).script3.voltage(ind_chg(end)) - ... 
              data(ind25).script3.voltage(ind_chg(end)+1); 
    
    % Bound resistant values into reasonable range
    % IR1_dis and IR2_chg are all under SOC=100%
    % IR2_dis and IR1_chg are all under SOC=0%
    % Here similar ones are used to bound each other
    IR1_dis = min(IR1_dis,2*IR2_chg); IR2_dis = min(IR2_dis,2*IR1_chg);
    IR1_chg = min(IR1_chg,2*IR2_dis); IR2_chg = min(IR2_chg,2*IR1_dis);
    
    %-----------------------------------------------------------------
    % Prepare for SOC-OCV relationship, intuition is create two arrays 
    % which contain SOC and OCV values correspondingly.
    %-----------------------------------------------------------------

    blend_dis = (0:length(ind_dis)-1)/(length(ind_dis)-1);
    % Assume a linear relationship between R and SOC, interpol R and SOC
    IRblend_dis = IR1_dis + (IR2_dis - IR1_dis)*blend_dis(:);
    % Compensate the voltage sharing from resistance to get OCV(disV)
    disV = data(ind25).script1.voltage(ind_dis) + IRblend_dis;
    disZ = 1 - data(ind25).script1.disAh(ind_dis)/Q_25;
    % Ensure the discharge SOC decrease from 1
    disZ = disZ + (1 - disZ(1));
    file_data(ind25).disZ = disZ;
    % store the terminal voltage rather the OCV
    file_data(ind25).disV = data(ind25).script1.voltage(ind_dis);

    blend_chg = (0:length(ind_chg)-1)/(length(ind_chg)-1);
    IRblend_chg = IR1_chg + (IR2_chg - IR1_chg)*blend_chg(:);
    chgV = data(ind25).script3.voltage(ind_chg) - IRblend_chg;
    chgZ = data(ind25).script3.chgAh(ind_chg)/Q_25;
    chgZ = chgZ - chgZ(1);
    file_data(ind25).chgZ = chgZ;
    file_data(ind25).chgV = data(ind25).script3.voltage(ind_chg);

    % compute OCV difference between charge and dischge at 50% SOC
    % force i*R compensated curve to pass half-way between each charge
    % and discharge at this point
    delta_V50 = interp1(chgZ, chgV, 0.5) - interp1(disZ, disV, 0.5);
    ind_chg50 = find(chgZ < 0.5);
    OCV_chg = chgV(ind_chg50) - delta_V50*chgZ(ind_chg50);
    Z_chg = chgZ(ind_chg50);

    ind_dis50 = find(disZ > 0.5);
    OCV_dis = flipud(disV(ind_dis50) + delta_V50*(1-disZ(ind_dis50)));
    Z_dis = flipud(disZ(ind_dis50));

    file_data(ind25).rawocv = interp1([Z_chg; Z_dis],[OCV_chg; OCV_dis],SOC,...
                               'linear','extrap');
    % "rawocv" now has our best guess of true OCV at this temperature
    file_data(ind25).temp = data(ind25).temp;
    
    



    % ------------------------------------------------------------------
    % Create SOC and OCV arrays under other temperatures
    % The only difference is coulombic efficiencies
    % ------------------------------------------------------------------
    for k = indnot25',
        data(k).script2.chgAh = data(k).script2.chgAh*eta_25;
        data(k).script4.chgAh = data(k).script4.chgAh*eta_25;
        eta(k) = (data(k).script1.disAh(end) + data(k).script2.disAh(end) + ...
                  data(k).script3.disAh(end) + data(k).script4.disAh(end) - ...
                  data(k).script2.chgAh(end) - data(k).script4.chgAh(end))/ ...
                  (data(k).script1.chgAh(end)+ data(k).script3.chgAh(end));
        data(k).script1.chgAh = eta(k)*data(k).script1.chgAh;         
        data(k).script3.chgAh = eta(k)*data(k).script3.chgAh;  

        tol_Q(k) = data(k).script1.disAh(end) + data(k).script2.disAh(end) - ...
            data(k).script1.chgAh(end) - data(k).script2.chgAh(end);

        ind_dis  = find(data(k).script1.step == 2);
        IR1_dis = data(k).script1.voltage(ind_dis(1)-1) - ... 
                  data(k).script1.voltage(ind_dis(1)); 
        IR2_dis = data(k).script1.voltage(ind_dis(end)+1) - ... 
                  data(k).script1.voltage(ind_dis(end)); 

        ind_chg = find(data(k).script3.step == 2);
        IR1_chg = data(k).script3.voltage(ind_chg(1)) - ... 
                  data(k).script3.voltage(ind_chg(1)-1); 
        IR2_chg = data(k).script3.voltage(ind_chg(end)) - ... 
                  data(k).script3.voltage(ind_chg(end)+1); 
        IR1_dis = min(IR1_dis,2*IR2_chg); IR2_dis = min(IR2_dis,2*IR1_chg);
        IR1_chg = min(IR1_chg,2*IR2_dis); IR2_chg = min(IR2_chg,2*IR1_dis);

        blend_dis = (0:length(ind_dis)-1)/(length(ind_dis)-1);
        IRblend_dis = IR1_dis + (IR2_dis - IR1_dis)*blend_dis(:);
        disV = data(k).script1.voltage(ind_dis) + IRblend_dis;
        disZ = 1 - data(k).script1.disAh(ind_dis)/Q_25;
        disZ = disZ + (1 - disZ(1));
        file_data(k).disZ = disZ;
        file_data(k).disV = data(k).script1.voltage(ind_dis);

        blend_chg = (0:length(ind_chg)-1)/(length(ind_chg)-1);
        IRblend_chg = IR1_chg + (IR2_chg - IR1_chg)*blend_chg(:);
        chgV = data(k).script3.voltage(ind_chg) - IRblend_chg;
        chgZ = data(k).script3.chgAh(ind_chg)/Q_25;
        chgZ = chgZ - chgZ(1);
        file_data(k).chgZ = chgZ;
        file_data(k).chgV = data(k).script3.voltage(ind_chg);
        
        delta_V50 = interp1(chgZ, chgV, 0.5) - interp1(disZ, disV, 0.5);
        ind_chg50 = find(chgZ < 0.5);
        OCV_chg = chgV(ind_chg50) - delta_V50*chgZ(ind_chg50);
        Z_chg = chgZ(ind_chg50);
        ind_dis50 = find(disZ > 0.5);
        OCV_dis = flipud(disV(ind_dis50) + delta_V50*(1-disZ(ind_dis50)));
        Z_dis = flipud(disZ(ind_dis50));
        file_data(k).rawocv = interp1([Z_chg; Z_dis],[OCV_chg; OCV_dis],SOC,...
                                   'linear','extrap');
        file_data(k).temp = data(k).temp;
    end
    model.OCVeta = eta;
    model.OCVQ = tol_Q;
    model.name = cellID;

    % Concatenate OCV and temperatures into two arrays
    OCV_tol = []; temp_tol = [];
    for i = 1:numtemps,
        if file_data(i).temp > 0,
            OCV_tol = [OCV_tol; file_data(i).rawocv];
            temp_tol = [temp_tol; file_data(i).temp];
        end
    end

    numtemps_kept = size(OCV_tol,1);

    % OCV(T,Z) = OCV0(Z) + T*OCVrel(T,Z)
    % Calculate the OCV-SOC(Z) relationship
    % Y = AX 
    % Y = [OCV(T1,Z_k),...,OCV(Tn,Z_k)]' = OCV_tol(:,k)
    % A = [1 T1;...;1 Tn]'; X = [OCV0(Z_k);OCVrel(Z_k)].
    OCV0 = zeros(size(SOC)); OCVrel = OCV0;
    A = [ones([numtemps_kept, 1]), temp_tol]; Y = OCV_tol;
    for i = 1:length(SOC),
        X = A\OCV_tol(:,i);
        OCV0(i) = X(1);
        OCVrel(i) = X(2);
    end
    model.OCV0 = OCV0;
    model.OCVrel = OCVrel;
    model.SOC = SOC;

    % ------------------------------------------------------------------
    % Make SOC0 and SOCrel
    % Do same kind of analysis to find SOC as a function of OCV
    % SOC(T,V) = SOC0(V) + T*SOCrel(T,V)
    % Y = AX 
    % Y = [SOC(T1,V_k),...,SOC(Tn,V_k)]' = SOC_tol(:,k)
    % A = [1 T1;...;1 Tn]'; X = [SOC0(V_k);SOCrel(V_k)].
    % ------------------------------------------------------------------
    Z = -0.1:0.01:1.1; % test soc vector
    V = minV-0.01:0.01:maxV+0.01;
    soc_tol = [];
    for T = filetemps',
    V1 = getOCVfromSOCTemp(Z,T,model); % V1 generated from OCV-SOC above
    soc_tol = [soc_tol; interp1(V1,Z,V)];
    end
    
    SOC0 = zeros(size(V)); SOCrel = SOC0; 
    A = [ones([numtemps,1]), filetemps]; 
    for k = 1:length(V),
    X = A\soc_tol(:,k); 
    SOC0(k) = X(1); 
    SOCrel(k) = X(2);
    end
    model.OCV = V;
    model.SOC0 = SOC0;
    model.SOCrel = SOCrel;
    

    %---------------------------------------------------------------------
    % Plot figures to show the OCV-SOC curves under different temperature
    % Plot individual OCV-SOC figure for 25 degC temperature with disV/chgV
    %---------------------------------------------------------------------
    figure;
    Legend = cell(length(filetemps),1);
    for i = 1:length(filetemps),
        plot(100*SOC, getOCVfromSOCTemp(SOC,file_data(i).temp, model))
        Legend{i} = strcat(int2str(file_data(i).temp),'degC');
        hold on
    end
    legend(Legend)
    xlabel('SOC')
    ylabel('Voltage')
    title('OCV-SOC under different temps')
    hold off
    
    figure;
    plot(100*SOC, getOCVfromSOCTemp(SOC,file_data(ind25).temp, model)*10,'r');
    hold on
    %plot(100*disZ, disV, 'b--', 'linewidth',1)
    %plot(100*chgZ, chgV, 'g--', 'linewidth',1)
    
    plot(100*file_data(ind25).disZ, file_data(ind25).disV*10, 'b--', 'linewidth',1);
    plot(100*file_data(ind25).chgZ, file_data(ind25).chgV*10, 'g--', 'linewidth',1);
    legend(["Estimated OCV", "Terminal disV", "Terminal chgV"])
    %legend(["Estimated OCV", "dis-ocv", "chg-ocv","Terminal disV", "Terminal chgV"])
    xlabel('SOC')
    ylabel('Voltage')
    title('OCV-SOC at 25 degC with disV/chgV')
    hold off

    






    



  

    

