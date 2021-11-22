% function generateDynamic
% The inputs: 
% - data: dynamic operation data for the selected battery
% - model: the OCV-SOC relationship gain from generateOCVSOC.m
% - num_rc: the number of R-C pairs
% The output:
% - model: A modified model, contains both static (OCV-SOC) and dynamic
% field model.


function model = generateDynamic(data, model, num_rc)
  % ------------------------------------------------------------------
  % Step 1: Compute capacity and coulombic efficiency for every test
  % ------------------------------------------------------------------
  tol_temps = [data(:).temp];
  tol_etas = 0*tol_temps;
  tol_Qs = 0*tol_temps;

  ind25 = find(tol_temps == 25);
  if isempty(ind25),
      error('Must have a test at 25degC');
  end
  indnot25 = find(tol_temps ~= 25);
  
  tol_disAh = data(ind25).script1.disAh(end) + ...
              data(ind25).script2.disAh(end) + ...
              data(ind25).script3.disAh(end);
  tol_chgAh = data(ind25).script1.chgAh(end) + ...
              data(ind25).script2.chgAh(end) + ...
              data(ind25).script3.chgAh(end);
  % 25 degC coulombic efficiency
  eta_25 = tol_disAh/tol_chgAh;
  data(ind25).eta = eta_25;
  tol_etas(ind25) = eta_25;
  data(ind25).script1.chgAh = data(ind25).script1.chgAh*eta_25;
  data(ind25).script2.chgAh = data(ind25).script2.chgAh*eta_25; 
  data(ind25).script3.chgAh = data(ind25).script3.chgAh*eta_25; 

  Q_25 = data(ind25).script1.disAh(end) + data(ind25).script2.disAh(end) - ...
         data(ind25).script1.chgAh(end) - data(ind25).script2.chgAh(end);
  data(ind25).Q = Q_25;
  tol_Qs(ind25) = Q_25;

  for k = indnot25,
    data(k).script2.chgAh = data(k).script2.chgAh*eta_25;
    data(k).script3.chgAh = data(k).script3.chgAh*eta_25;
    eta = (data(k).script1.disAh(end) + data(k).script2.disAh(end)+...
           data(k).script3.disAh(end) - data(k).script2.chgAh(end)-...
           data(k).script3.chgAh(end))/data(k).script1.chgAh(end);
    data(k).script1.chgAh = eta*data(k).script1.chgAh;
    data(k).eta = eta; tol_etas(k) = eta;
    
    Q = data(k).script1.disAh(end) + data(k).script2.disAh(end) - ...
          data(k).script1.chgAh(end) - data(k).script2.chgAh(end);
    data(k).Q = Q; tol_Qs(k) = Q;
  end
    
  % there may exist different data sets from same testing temperature
  % and the toltal capacity Q and coulumbic efficiency eta should be same
  % under same temperature.
  model.temps = unique(tol_temps); num_temps = length(model.temps);
  
  model.etaParam = NaN(1, num_temps);
  model.QParam = NaN(1, num_temps);
  for i = 1:num_temps,
      model.etaParam(i) = mean(tol_etas(tol_temps == model.temps(i)));
      model.QParam(i) = mean(tol_Qs(tol_temps == model.temps(i)));
  end

  % ------------------------------------------------------------------
  % Step 2: Compensate OCV in charging parts
  % ------------------------------------------------------------------

  for i = 1:length(data),
      curr_eta = model.etaParam(i);
      eta_I = data(i).script1.current;
      % I < 0 means charging, and change capacity unit (Ah->As*3600)
      eta_I(eta_I < 0) = eta_I(eta_I < 0)*curr_eta;
      data(i).Z = 1 - cumsum([0, eta_I(1:end-1)])/(data(i).Q*3600);
      data(i).OCV = getOCVfromSOCTemp(data(i).Z, tol_temps(i), model);
  end

  % Now begin to optimize for dynamic parameters
  model.R0Param = NaN(1,num_temps); % "R0" ohmic resistance parameter
  model.RCParam = NaN(num_temps,num_rc); % time const.
  model.RParam  = NaN(num_temps,num_rc); % Rk

  for ind_temp = 1:num_temps,
      model.GParam(ind_temp) = 0;
      curr_temp = model.temps(ind_temp);
      fprintf('Processing temperature %d\n', curr_temp);
      
      curr_ind = find(tol_temps == curr_temp);
      num_files = length(curr_ind);

      xplots = ceil(sqrt(num_files));
      yplots = ceil(num_files/xplots);
      rmserr = zeros(1,xplots*yplots);
      Q = abs(getParamESC('QParam',curr_temp,model));
      eta = abs(getParamESC('etaParam',curr_temp,model));
      RC = getParamESC('RCParam',curr_temp,model);
      num_rc = length(RC);

      for curr_file = 1:num_files,
        ik = data(curr_ind(curr_file)).script1.current(:);
        vk = data(curr_ind(curr_file)).script1.voltage(:);
        tk = (1:length(vk))-1;
        etaik = ik; etaik(ik<0) = etaik(ik<0)*eta;
        
        %h=0*ik; sik = 0*ik;
        %fac=exp(-abs(0*etaik/(3600*Q)));
        %for k=2:length(ik),
        %    h(k)=fac(k-1)*h(k-1)-(1-fac(k-1))*sign(ik(k-1));
        %    sik(k) = sign(ik(k));
        %    if abs(ik(k))<Q/100, sik(k) = sik(k-1); end
        %end

        % First modeling step: Compute error with model = OCV only
        
        vest1 = data(curr_ind(curr_file)).OCV;
        
        verr = vk - vest1;

        % Second modeling step: Compute time constants in "A" matrix
        % This technique is called subspace system identification which is
        % used to replace nonlinear optimization. As it can finds solution
        % in one step.
        np = num_rc;
        while 1,
            A = SISOsubid(-diff(verr),diff(etaik),np);
            eigA = eig(A); 
            eigA = eigA(eigA == conj(eigA));  % make sure real
            eigA = eigA(eigA > 0 & eigA < 1); % make sure in range
            okpoles = length(eigA); np = np+1;
            if okpoles >= num_rc, break; end
            fprintf('Trying np = %d\n',np);
        end
        RCfact = sort(eigA); RCfact = RCfact(end-num_rc+1:end);
        RC = -1./log(RCfact);
        % Simulate the R-C filters to find R-C currents
        % dlsim.m is in the control-system toolbox
        if exist('dlsim.m', 'file')
            vrcRaw = dlsim(diag(RCfact),1-RCfact,...
                     eye(num_rc),zeros(num_rc,1),etaik);
        else
            % a slower workaround if no control-system toolbox
            vrcRaw = zeros(length(RCfact),length(etaik));
            for vrcK = 1:length(etaik)-1,
            vrcRaw(:,vrcK+1) = diag(RCfact)*vrcRaw(:,vrcK)+(1-RCfact)*etaik(vrcK);
            end
            vrcRaw = vrcRaw';
        end

        H = [-etaik,-vrcRaw]; 
        W = H\verr;    
        R0 = W(1); Rfact = W(2:end)';
        ind = find(model.temps == data(curr_ind(curr_file)).temp,1);
        model.R0Param(ind) = R0;
        model.RCParam(ind,:) = RC';
        model.RParam(ind,:) = Rfact';

        vest2 = vest1 - R0*etaik - vrcRaw*Rfact';
        verr = vk - vest2;

        % plot voltages for 25 degC
        if curr_temp == 25,
            figure(1); 
            subplot(yplots,xplots,curr_file); 
            %plot(tk(1:10:end)/60,vk(1:10:end),tk(1:10:end)/60,...
            % vest1(1:10:end),tk(1:10:end)/60,vest2(1:10:end));  
            plot(tk(1:10:end)/60,vk(1:10:end),'b');
            hold on
            plot(tk(1:10:end)/60,vest2(1:10:end),'r');
            plot(tk(1:10:end)/60,vest1(1:10:end),'g');
            hold off
            xlabel('Time (min)'); ylabel('Voltage (V)'); 
            title(sprintf('Voltage and estimates at T=%d',...
                      data(ind(curr_file)).temp));
            legend('voltage','vest2 (DYN)','vest1 (OCV)','location','southwest');
            
            % plot modeling errors
            figure(2); subplot(yplots,xplots,curr_file); 
            thetitle=sprintf('Modeling error at T = %d',data(ind(curr_file)).temp);
            plot(tk(1:10:end)/60,verr(1:10:end)); title(thetitle);
            xlabel('Time (min)'); ylabel('Error (V)');
            ylim([-0.1 0.1]); 
            drawnow
        end

        % Compute RMS error only on data roughly in 5% to 95% SOC
        v1 = getOCVfromSOCTemp(0.95,data(ind(curr_file)).temp,model);
        v2 = getOCVfromSOCTemp(0.05,data(ind(curr_file)).temp,model);
        N1 = find(vk<v1,1,'first'); N2 = find(vk<v2,1,'first');
        if isempty(N1), N1=1; end; if isempty(N2), N2=length(verr); end
        rmserr(curr_file)=sqrt(mean(verr(N1:N2).^2));
      end
        
      cost = sum(rmserr);
      fprintf('RMS error = %0.2f (mV)\n',cost*1000);
      if isnan(cost), stop, end
      
      if curr_temp == 25,
        figure(3); theXlim = [min(model.temps) max(model.temps)];
        subplot(2,2,1); 
        plot(model.temps,model.QParam); 
        title('Capacity (Ah)'); 
        xlim(theXlim);
        subplot(2,2,2); 
        plot(model.temps,1000*model.R0Param); 
        title('Resistance (m\Omega)');
        subplot(2,2,3); 
        plot(model.temps,getParamESC('RCParam', model.temps,model));
        title('RC Time Constant (tau)');
        xlim(theXlim);
        subplot(2,2,4); 
        plot(model.temps,1000*getParamESC('RParam',model.temps,model));
        title('R in RC (m\Omega)');
        xlim(theXlim);
      end
  end


  

























  

  
   