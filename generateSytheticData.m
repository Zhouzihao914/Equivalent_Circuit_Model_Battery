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


tol_data_num = 100;
repeat_times = 50;
num_rc = 1;
temp = 25;
do_Qtime = 1;
dis_unit = udds_current/(8);
Qmax = model.QParam(find(model.temps == 25));
Q0 = Qmax;

% set up the hyper-params for a,b
u_a = 0.01/72000; sigma_a = 0.01/72000;
u_b = 0.01/3600; sigma_b = 0.01/3600;

% set up the cycle unit profile
[v,rc,z_dis,OCV,Qend] = ECMcell(dis_unit,temp,1,...
                            model,1,zeros(num_rc,1),Q0,u_a,1,u_b);

delta_q = (1-z_dis(end))*model.QParam(find(model.temps == temp));
chg_time = ceil(delta_q*3600/1);
chg_unit = -1*ones(chg_time,1);
short_rest_unit = zeros([300,1]);
long_rest_unit = zeros([1000,1]);
%all_curr_unit = [dis_unit;short_rest_unit;chg_unit;long_rest_unit];
%all_curr_unit = all_curr_unit(1:7200);

for i = 1:tol_data_num,
    % ensure the a,b value is reasonable
    curr_a = max(0.1*u_a, normrnd(u_a, sigma_a));
    curr_a = min(curr_a,10*u_a);

    curr_b = max(0.1*u_b,normrnd(u_b, sigma_b));
    curr_b = min(curr_b,1.5*u_b);

    all_data(i).a = curr_a;
    all_data(i).b = curr_b;

    Q0 = Qmax;
    z0 = 1;

    cycle_v = [];
    cycle_curr = [];
    cycle_z = [];
    cycle_Q = [];
    cycle_time = 0;
    Q_q = [];
    Q_time = [];
    real_chg_time = [];
    
    
    Q_r = model.R0Param(find(model.temps == 25)) * ones(repeat_times,1);   
    
    all_data(i).Qinit = Q0;
    for ii = 1:repeat_times,
        cycle_num = ii;
        %discharge stage
        [vdis,rc,zdis,OCV,Qdis] = ECMcell(dis_unit,temp,1,...
                            model,z0,zeros(num_rc,1),Q0,curr_a,1,curr_b);
        cycle_v = [cycle_v;vdis];
        cycle_curr = [cycle_curr;dis_unit];
        cycle_z = [cycle_z;zdis];
        cycle_Q = [cycle_Q;Qdis];

        %short rest stage
        [vrest1,rc,zrest1,OCV,Qrest1] = ECMcell(short_rest_unit,temp,1,...
                            model,zdis(end),zeros(num_rc,1),Qdis(end),curr_a,0,curr_b);
        cycle_v = [cycle_v;vrest1];
        cycle_curr = [cycle_curr;short_rest_unit];
        cycle_z = [cycle_z;zrest1];
        cycle_Q = [cycle_Q;Qrest1];

        %charge stage
        [vchg,rc,zchg,OCV,Qchg] = ECMcell(chg_unit,temp,1,model,...
                        zrest1(end),zeros(num_rc,1),Qrest1(end),curr_a,0,curr_b);

        %check if the voltage beyond cutoff volt
        vcut = find(vchg >= 3.6);
        if length(vcut) > 0,
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
        
        %long rest stage
        [vrest2,rc,zrest2,OCV,Qrest2] = ECMcell(long_rest_unit,temp,1,model,...
                              zchg(end),zeros(num_rc,1),Qchg(end),curr_a,0,curr_b);

        Q0 = Qrest2(end);
        cycle_v = [cycle_v; vrest2];
        cycle_curr = [cycle_curr; long_rest_unit];
        cycle_z = [cycle_z; zrest2];
        cycle_Q = [cycle_Q; Qrest2];
        cycle_time = cycle_time + length(cycle_v);

        Q_q = [Q_q; Qrest2(end)];
        Q_time = [Q_time; cycle_time];
    end

    Q_cycle = [1:1:repeat_times]';
    cycle_time = [1:1:length(cycle_v)]';
    cycle_temp = temp*ones(length(cycle_v),1);
    cycle_table = [cycle_time, cycle_curr, cycle_v, cycle_temp];
    
    q_table = [Q_time, Q_cycle, Q_q, Q_r];
    
    all_q{i} = array2table(q_table,'VariableNames',q_table_variables);
    all_data(i).cycles = array2table(cycle_table,'VariableNames',cycle_table_variables);
    all_data(i).soc = cycle_z;
    all_data(i).minsoc = min(cycle_z);
    all_data(i).maxsoc = max(cycle_z);
    all_data(i).Qend = cycle_Q;
    all_data(i).chgtime = real_chg_time';
end

Q = all_q';
cycledata = all_data';

save('synthetic_Q.mat', "Q", '-v7.3');
save('synthetic_cycledata.mat', "cycledata", '-v7.3');

figure(1)
for j = 1:1:tol_data_num,
    t = [1:1:length(all_data(j).cycles.volt)];
    plot(t/3600, all_data(j).cycles.volt);
    hold on
end
xlabel('Time(hr)')
ylabel('Voltage(volt)')
title('Voltage Curve')
hold off

figure(2)
for j = 1:1:tol_data_num,
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
    t = [1:1:length(all_data(j).cycles.volt)];
    plot(t/3600, all_data(j).soc);
    hold on
end
xlabel('Time(hr)')
ylabel('SOC(%)')
title('SOC Curve')
hold off

figure(4)
s = {};
for j = 1:1:tol_data_num,
    t = [1:1:length(all_data(j).cycles.volt)];
    plot(t/3600, all_data(j).Qend);
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
for j = 1:1:tol_data_num,
    c = [1:1:repeat_times];
    plot(c, all_data(j).chgtime);
    hold on
end
xlabel('Cycles')
ylabel('Charging Time(hr)')
title('Charge Time')
hold off

figure(6)
for j = 1:1:tol_data_num,
    t = [1:1:length(all_data(j).cycles.volt)];
    plot(t/3600, all_data(j).Qend .* all_data(j).soc);
    hold on
end
xlabel('Time(hr)')
ylabel('Coulomb Quantity(Ah)')
title('Real Coulomb Quantity')
hold off


