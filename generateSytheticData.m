clear
close all
load 'udds/udds_curr.mat'

modelFile = 'A123model-ocv.mat';
load(modelFile);

format long
%chg_curr_value = (1-0.0253)*2.0495/0.9944/2;
repeat_times = 50;
num_rc = 1;
temp = 25;
do_Qtime = 1;
dis_unit = udds_current/10;
Qmax = model.QParam(find(model.temps == 25));
Q0 = Qmax;

[v,rc,z_dis,OCV,Qend] = ECMcell(dis_unit,temp,1,...
                            model,1,zeros(num_rc,1),Q0,0);

delta_q = (1-z_dis(end))*model.QParam(find(model.temps == temp));
chg_time = floor(delta_q*3600/1);
chg_unit = -ones([chg_time,1]);
short_rest_unit = zeros([300,1]);
long_rest_unit = zeros([3981,1]);

all_unit = [dis_unit;short_rest_unit;chg_unit;long_rest_unit];
z0 = 1;
tol_v = [];
tol_curr = [];
tol_z = [];

for i = 1:repeat_times,
    cycle_num = i;
    [v,rc,z,OCV,Qend] = ECMcell(all_unit,temp,1,...
                            model,z0,zeros(num_rc,1),Q0,cycle_num);
    Q0 = Qend;
    tol_v = [tol_v;v];
    tol_curr = [tol_curr;all_unit];
    tol_z = [tol_z;z];
end

figure(1);
plot([1:1:length(tol_curr)]/3600, tol_curr);
xlabel('Time(hr)')
ylabel('Current(A)')
%ylim([0,1])

figure(2);
plot([1:1:length(tol_curr)]/3600, tol_v);
xlabel('Time(hr)')
ylabel('Voltage(V)')

figure(3);
plot([1:1:length(tol_curr)]/3600, tol_z);
xlabel('Time(hr)')
ylabel('SOC(%)')
%ylim([0,1])

