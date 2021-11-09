% script onesampleOCV.m
% load the data from a cell <A123> tested under eight different 
% temperatures. This script calls getOCV.m to create OCV-SOC
% relationsip, then saves the model to a model file.

% This script is a self-reproduced simplified version from the one provided
% by Prof. Plett, Gregory L. <http://mocha-java.uccs.edu/BMS1/index.html>

clear all; close all;
cellID = 'A123';
temps = [-25 -15 -5 5 15 25 35 45];
minV = 2.0; maxV = 3.75;
OCVDir = sprintf('%s_OCV',cellID);
data = zeros([0 length(temps)]);

%--------------------------------------------------------------------
% load cell test data 
%--------------------------------------------------------------------
for k = 1:length(temps),
    if temps(k) < 0,    
        filename = sprintf('%s/%s_OCV_N%02d.mat', ...
            OCVDir,cellID,abs(temps(k)));
    else
        filename = sprintf('%s/%s_OCV_P%02d.mat', ... 
            OCVDir,cellID,temps(k));
    end
    load(filename);
    data(k).temp = temps(k);
    data(k).script1 = OCVData.script1;
    data(k).script2 = OCVData.script2;
    data(k).script3 = OCVData.script3;
    data(k).script4 = OCVData.script4;
    
    %plot the voltage curve under temperature 25
    if data(k).temp == 25,
        figure
        allTimes = [data(k).script1.time; data(k).script2.time;
                 data(k).script3.time; data(k).script4.time];
        allVols = [data(k).script1.voltage; data(k).script2.voltage;
                 data(k).script3.voltage; data(k).script4.voltage];
        script_time1 = data(k).script1.time;
        script_time2 = data(k).script1.time(end) + data(k).script2.time;
        script_time3 = script_time2(end) + data(k).script3.time;
        script_time4 = script_time3(end) + data(k).script4.time;
        plot(script_time1/3600, data(k).script1.voltage,'b')
        hold on
        plot(script_time2/3600, data(k).script2.voltage,'r');
        plot(script_time3/3600, data(k).script3.voltage,'b');
        plot(script_time4/3600, data(k).script4.voltage,'r');
        xlabel('Time(hr)')
        ylabel('Voltage(volt)')
        title('Voltage Curves at 25 degC');
        hold off
    end    
end

% then, call "generateOCVSOC" to do the actual data processing
model = generateOCVSOC(data,cellID,minV,maxV);
save(sprintf('%smodel-ocv.mat',cellID),'model'); % save model file

