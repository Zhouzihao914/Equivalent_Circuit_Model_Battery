% Simulates synthetic battery usage data 
% The UDDS(Urban Dynamometer Driving Schedule) is used as dynamic profile
% for battery usage. Also the battery discharge capacity drop due to
% battery aging is included.
%
% 'func_id' is used to customize aging function for discharge capacity.
% (0-no aging effect; 1-linear aging; 2-exponential aging)
% 
% Running generateOCVSOC.m and generateDynamic.m
% before running this script

clear
close all
load 'udds/udds_curr.mat'

modelFile = 'A123model-ocv.mat';
load(modelFile);

format long
all_data = struct;
all_q = {};
cycle_table_variables = {'time','curr','volt','temp'};
q_table_variables = {'time','cycles','capacity','resistance'};

% Expected sample size and repeated usage cycles for each sample. It may
% takes more than 10 minutes if these parameters are large (~100/~50)
% also the avaliable repeat_times is closely related to your aging func, as
% the battery may aging to end of life before your predefined repeat_times.
tol_data_num = 5;
repeat_times = 60;
func_id = 2;

% ECM parameters
num_rc = 1;
temp = 25;
do_Qtime = 1;
dis_unit = udds_current/(8);
Qmax = model.QParam(find(model.temps == 25));
Q0 = Qmax;

% set up the hyper-params for a,b (linear disQ aging)
% u_a = 0.001/72000; sigma_a = 0.001/72000;
% u_b = 0.01/3600; sigma_b = 0.01/3600;

% set up the hyper-params for a,b (Exp disQ aging)
u_a = 0.001/54000; sigma_a = 0.001/54000;
u_b = 6; sigma_b = 3;

% set up the cycle unit profile
[v,rc,z_dis,OCV,Qend] = ECMcell(dis_unit,temp,1,...
                            model,1,zeros(num_rc,1),Q0,0,u_a,u_b,0);
delta_q = (1-z_dis(end))*model.QParam(find(model.temps == temp));
chg_time = ceil(delta_q/1)*3600;
chg_unit = -1*ones(chg_time,1);
short_rest_unit = zeros([300,1]);
long_rest_unit = zeros([1000,1]);

for i = 1:tol_data_num
    % ensure the a,b value is reasonable
    if i == 1
        curr_a = u_a;
        curr_b = u_b;
    else
        curr_a = max(0.1*u_a, normrnd(u_a, sigma_a));
        curr_a = min(curr_a,5*u_a);
        curr_b = max(0.1*u_b,normrnd(u_b, sigma_b));
        curr_b = min(curr_b,2*u_b);
    end
    
    %for each sample, cell is initialized as 100% SOC and Qmax (new cell)
    all_data(i).a = curr_a;
    all_data(i).b = curr_b;
    Q0 = Qmax;
    z0 = 1;
    
    % store simulation data for later usage
    cycle_v = [];
    cycle_curr = [];
    cycle_z = [];
    cycle_Q = [];
    cycle_time = 0;
    cycle_index = [];
    real_chg_time = [];
    
    all_data(i).Qinit = Q0;
    for ii = 1:repeat_times,
        cycle_num = ii;
        %discharge stage
        [vdis,~,zdis,~,Qdis] = ECMcell(dis_unit,temp,1,...
                            model,z0,zeros(num_rc,1),Q0,func_id,curr_a,curr_b,ii);

        cycle_v = [cycle_v;vdis];
        cycle_curr = [cycle_curr;dis_unit];
        cycle_z = [cycle_z;zdis];
        cycle_Q = [cycle_Q;Qdis];
        cycle_index = [cycle_index; ii*ones([length(vdis),1])];

        %short rest stage
        [vrest1,~,zrest1,~,Qrest1] = ECMcell(short_rest_unit,temp,1,...
                            model,zdis(end),zeros(num_rc,1),Qdis(end),func_id,curr_a,0,ii);
        cycle_v = [cycle_v;vrest1];
        cycle_curr = [cycle_curr;short_rest_unit];
        cycle_z = [cycle_z;zrest1];
        cycle_Q = [cycle_Q;Qrest1];
        cycle_index = [cycle_index; ii*ones([length(vrest1),1])];

        %charge stage
        [vchg,~,zchg,~,Qchg] = ECMcell(chg_unit,temp,1,model,...
                        zrest1(end),zeros(num_rc,1),Qrest1(end),func_id,curr_a,0,ii);

        %check if the voltage beyond cutoff volt
        vcut = find(vchg >= 3.6);
        if ~isempty(vcut)
            vchg = vchg(1:vcut(1)-1);
            zchg = zchg(1:vcut(1)-1);
            Qchg = Qchg(1:vcut(1)-1);
            curr_chg = chg_unit(1:vcut(1)-1);
        else
            curr_chg = chg_unit;

        end
        real_chg_time = [real_chg_time;length(curr_chg)];
        cycle_v = [cycle_v;vchg];
        cycle_curr = [cycle_curr;curr_chg];
        cycle_z = [cycle_z;zchg];
        cycle_Q = [cycle_Q;Qchg];
        cycle_index = [cycle_index; ii*ones([length(vchg),1])];

        %long rest stage
        [vrest2,rc,zrest2,OCV,Qrest2] = ECMcell(long_rest_unit,temp,1,model,...
                              zchg(end),zeros(num_rc,1),Qchg(end),func_id,curr_a,0,ii);

        Q0 = Qrest2(end);
        cycle_v = [cycle_v; vrest2];
        cycle_curr = [cycle_curr; long_rest_unit];
        cycle_z = [cycle_z; zrest2];
        cycle_Q = [cycle_Q; Qrest2];
        cycle_time = cycle_time + length(cycle_v);

        %Q_q = [Q_q; Qrest2];
        %Q_time = [Q_time; cycle_time];
        cycle_index = [cycle_index; ii*ones([length(vrest2),1])];
    end
    
    cycle_time = [1:1:length(cycle_v)]';
    cycle_temp = temp*ones(length(cycle_v),1);
    cycle_table = [cycle_time, cycle_curr, cycle_v, cycle_temp];
    
    Q_r = model.R0Param(find(model.temps == 25)) * ones(length(cycle_v),1);
    q_table = [cycle_time, cycle_index, cycle_Q, Q_r];
    
    all_q{i} = array2table(q_table,'VariableNames',q_table_variables);
    all_data(i).cycles = array2table(cycle_table,'VariableNames',cycle_table_variables);
    all_data(i).soc = cycle_z;
    all_data(i).minsoc = min(cycle_z);
    all_data(i).maxsoc = max(cycle_z);
    all_data(i).Q = cycle_Q;
    all_data(i).chgtime = real_chg_time';
end

Q = all_q';
cycledata = all_data';

save('synthetic_Q_exp.mat', "Q", '-v7.3');
save('synthetic_cycledata_exp.mat', "cycledata", '-v7.3');


figure(1)
for j = 1:1:tol_data_num
    t = [1:1:length(all_data(j).cycles.volt)];
    plot(t/3600, all_data(j).cycles.volt);
    hold on
end
xlabel('Time(hr)')
ylabel('Voltage(volt)')
title('Voltage Curve')
hold off

figure(2)
for j = 1:1:tol_data_num
    t = [1:1:length(all_data(j).cycles.volt)];
    plot(t/3600, all_data(j).cycles.curr);
    hold on
end
xlabel('Time(hr)')
ylabel('Current(A)')
title('Current Curve')
hold off

figure(3)
for j = 1:1:tol_data_num,
    plot(all_data(j).cycles.time/3600, all_data(j).soc);
    hold on
end
xlabel('Time(hr)')
ylabel('SOC(%)')
title('SOC Curve')
hold off

figure(4)
s = {};
for j = 1:1:tol_data_num
    plot(all_data(j).cycles.time/3600, all_data(j).Q);
    s{j} =  strcat('a: ', num2str(all_data(j).a,6),...
            ' b: ', num2str(all_data(j).b,4));
    hold on
end
legend(s)
xlabel('Time(hr)')
ylabel('disQ(Ah)')
title('Discharge Capacity')
hold off

figure(5)
for j = 1:1:tol_data_num
    plot(all_data(j).cycles.time/3600, all_data(j).Q .* all_data(j).soc);
    s{j} =  strcat('a: ', num2str(all_data(j).a,6),...
            ' b: ', num2str(all_data(j).b,4));
    hold on
end
legend(s)
xlabel('Time(hr)')
ylabel('Coulomb Quantity(Ah)')
title('Real Coulomb Quantity')
hold off


