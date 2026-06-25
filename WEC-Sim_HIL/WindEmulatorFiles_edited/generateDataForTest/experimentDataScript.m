clearvars; close all; clc;

wMin_rpm = -200;%[rpm] Minimum shaft speed to be tested
wMax_rpm = -1600;%[rpm] Maximum shaft speed to be tested
wIncrement_rpm = -200;%[rpm] Increments for shaft speed tests
wShaft = linspace(wMin_rpm,wMax_rpm,(wMax_rpm-wMin_rpm)/wIncrement_rpm+1);

torqueMin_Nm = -10;%[Nm] Minimum torque to be tested
torqueMax_Nm = -70;%[Nm] Maximum torque to be tested
torqueIncrement_Nm = -10;% Increments for torque tests
torqueShaft = linspace(torqueMin_Nm,torqueMax_Nm,(torqueMax_Nm-torqueMin_Nm)/torqueIncrement_Nm+1);

freq = 250;%Data frequency
tTest = 40;%[s] Time for each test

% Create the .mat file
experimentData = experimentDataFunc(wShaft,torqueShaft,freq,tTest);

% time vector used for timeseries object
time = 0:1/freq:(1/freq)*(length(experimentData.wShaft_rpm)-1);

experimentData.time = time;
% Uncomment the following line to save the data in a .mat file

save('experimentDataB.mat','experimentData')


figure(1)
subplot(3,1,1)
plot(time, experimentData.wShaft_rpm)
ylabel('Shaft Speed [rpm]')
xlim([0 max(time)])

subplot(3,1,2)
plot(time, experimentData.torqueShaft_Nm)
ylabel('Torque [Nm]')
xlim([0 max(time)])

subplot(3,1,3)
plot(experimentData.count)
ylabel('count')
xlim([0 max(time)])