
load 'current profile unit A.mat'

modelFile = 'A123model-ocv.mat';
load(modelFile);

repeat_times = 30;
sythetic_currents = [];
for i = 1:repeat_times,
sythetic_currents = [sythetic_currents;curr_profile_unit'];
end

num_rc = 1;
do_Qtime = 1;
temp = 25;
[vk,rck,zk,OCV] = ECMcell(sythetic_currents,temp,1,...
                            model,1,zeros(num_rc,1),do_Qtime);

figure;
plot([1:1:length(sythetic_currents)]/3600, vk);

figure;
plot([1:1:length(sythetic_currents)]/3600, zk);
