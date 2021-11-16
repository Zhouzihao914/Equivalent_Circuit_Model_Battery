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
repeat_times = 100;
num_rc = 1;
temp = 25;
do_Qtime = 1;
dis_unit = udds_current/(8);
Qmax = model.QParam(find(model.temps == 25));
Q0 = Qmax;

% set up the hyper-params for a,b
u_a = 0.0001; sigma_a = 0.00005;
u_b = 0.001; sigma_b = 0.0005;

% set up the cycle unit profile
[v,rc,z_dis,OCV,Qend] = ECMcell(dis_unit,temp,1,...
                            model,1,zeros(num_rc,1),Q0,u_a,u_b);

delta_q = (1-z_dis(end))*model.QParam(find(model.temps == temp));
chg_time = floor(delta_q*3600/1);
chg_unit = -ones([chg_time,1]);
short_rest_unit = zeros([300,1]);
long_rest_unit = zeros([4000,1]);
all_curr_unit = [dis_unit;short_rest_unit;chg_unit;long_rest_unit];
all_curr_unit = all_curr_unit(1:7200);

for i = 1:tol_data_num,
    % ensure the a,b value is reasonable
    curr_a = max(0.00001, normrnd(u_a, sigma_a));
    curr_a = min(curr_a,0.001);

    curr_b = max(0.0001,normrnd(u_b, sigma_b));
    curr_b = min(curr_b,0.01);

    all_data(i).a = curr_a;
    all_data(i).b = curr_b;

    Q0 = Qmax;
    z0 = 1;

    cycle_v = [];
    cycle_curr = [];
    cycle_z = [];

    Q_q = [];
    Q_r = model.R0Param(find(model.temps == 25)) * ones(repeat_times,1);   


    all_data(i).Qinit = Q0;
    for ii = 1:repeat_times,
        cycle_num = ii;
        [v,rc,z,OCV,Qend] = ECMcell(all_curr_unit,temp,1,...
                            model,z0,zeros(num_rc,1),Q0,curr_a,curr_b);
        Q0 = Qend;
        
        cycle_v = [cycle_v;v];
        cycle_curr = [cycle_curr;all_curr_unit];
        cycle_z = [cycle_z;z];
        Q_q = [Q_q;Qend];
    end
    
    Q_time = [1:1:repeat_times]'*length(all_curr_unit);
    Q_cycle = [1:1:repeat_times]';

    cycle_time = [1:1:repeat_times*length(all_curr_unit)]';
    cycle_temp = temp*ones(repeat_times*length(all_curr_unit),1);
    cycle_table = [cycle_time, cycle_curr, cycle_v, cycle_temp];
    
    q_table = [Q_time, Q_cycle, Q_q, Q_r];
    
    all_q{i} = array2table(q_table,'VariableNames',q_table_variables);
    all_data(i).cycles = array2table(cycle_table,'VariableNames',cycle_table_variables);
    all_data(i).soc = cycle_z;
    all_data(i).minsoc = min(cycle_z);
    all_data(i).maxsoc = max(cycle_z);
    all_data(i).Qend = Qend;
end

Q = all_q';
cycledata = all_data';
save('synthetic_Q.mat', "Q", '-v7.3');
save('synthetic_cycledata.mat', "cycledata", '-v7.3');
