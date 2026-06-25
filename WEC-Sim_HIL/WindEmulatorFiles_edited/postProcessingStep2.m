% post processing script to retrieve data from the target

clearvars; close all; clc;

%% === parameters =========================================================
%dataFile = 'windEmulatorStep4_2022_12_05_13_08_50';

%dataFile = 'windEmulatorStep4_2022_12_06_11_25_55'; % shake down test; confirmed all ABB side scale factors
%dataFile = 'windEmulatorStep4_2022_12_06_10_23_01'; % first partial SID matrix; stopped, as we hit torque limit at 90% (72Nm)
%dataFile = 'windEmulatorStep4_2022_12_06_10_46_29'; % first partial SID matrix with 97.5% torque limit; ran 200, 400, 600 rpm with 10...70Nm

%dataFile = 'windEmulatorStep4_2022_12_06_17_19_04.mat'; % short run to confirm file writing
%dataFile = 'windEmulatorStep4_2022_12_06_18_16_59.mat'; % 1st 15 cases of SID matrix

%dataFile = 'windEmulatorStep4_2022_12_07_09_42_49.mat';

%dataFile = 'windEmulatorStep4_2022_12_07_11_38_46.mat'; % re-tuning PI controller after changing filter times from ~500ms to 10ms (torque) and 3ms (speed) on both ABB drives
%dataFile = 'windEmulatorStep4_2022_12_07_11_57_59.mat'; % full Q1 SID matrix with updated drive filter settings and encoder-based control
%dataFile = 'windEmulatorStep4_2022_12_07_12_59_46.mat'; % full Q3 SID matrix; as Q1 above

%dataFile = 'windEmulatorStep4_2022_12_09_10_14_37.mat'; % full Q3 SID matrix; with HBM
dataFile = 'windEmulatorStep4_2024_01_09_14_44_10.mat'; % full Q3 SID matrix; with HBM


baseDir = fullfile('C:','simulink_processed');
dataFile = fullfile(baseDir,dataFile);

% =========================================================================

data = load(dataFile);

eCATBusName = 'eCATStatus';

for n = 1:data.data.numElements 

    if contains(string(data.data{n}.BlockPath.convertToCell),'hptoCtrl')
        indexHptoCtrl = n;
    end
    
    if contains(string(data.data{n}.BlockPath.convertToCell),'hptoSignals')
        indexHptoSignals = n;
    end
    
    if contains(string(data.data{n}.BlockPath.convertToCell),'acs880Signals')
        indexAcs880Signals = n;
    end      
    
    if contains(string(data.data{n}.BlockPath.convertToCell),'expCtrlSignals')
        indexExpCtrlSignals = n;
    end      
    
    if contains(string(data.data{n}.BlockPath.convertToCell),'invPowerAcs800')
        indexInvPowerAcs800Signals = n;
    end      
    
    if contains(string(data.data{n}.BlockPath.convertToCell),'invPowerAcs880')
        indexInvPowerAcs880Signals = n;
    end  
    
    if contains(string(data.data{n}.BlockPath.convertToCell),'shaftSignals')
        indexShaftSignals = n;
    end   
    
    if contains(string(data.data{n}.BlockPath.convertToCell),'acs800Signals')
        indexAcs800Signals = n;
    end   
    
    if contains(string(data.data{n}.BlockPath.convertToCell),'sidInfoSignals')
        indexSidInfoSignals = n;
    end   
    
end

indexInvPowerAcs880Signals = 9;

if ~exist('indexExpCtrlSignals','var')
    error('Did not find exp control bus variable in results file');
end


if ~exist('indexHptoCtrl','var')
    error('Did not find HPTO control bus variable in results file');
end

if ~exist('indexHptoSignals','var')
    error('Did not find MotionCtrl bus variable in results file');
end

if ~exist('indexSidInfoSignals','var')
    error('Did not find SID info signals variable in results file');
end


expCtrlData = data.data{indexExpCtrlSignals};
hptoCtrlData = data.data{indexHptoCtrl};
hptoSignalsData = data.data{indexHptoSignals};
invPowerAcs800Data = data.data{indexInvPowerAcs800Signals};
invPowerAcs880Data = data.data{indexInvPowerAcs880Signals};
shaftSignalsData = data.data{indexShaftSignals};
acs880SignalsData = data.data{indexAcs800Signals};
sidInfoSignalsData = data.data{indexSidInfoSignals};
% =========================================================================


%% === extract relevant signals ===========================================
time = expCtrlData.Values.ramp.Time;
expTime = squeeze(expCtrlData.Values.time.Data);
runCounter = squeeze(expCtrlData.Values.runCounter.Data);

acs800L1V = squeeze(invPowerAcs800Data.Values.L1Voltage.Data);
acs800L2V = squeeze(invPowerAcs800Data.Values.L2Voltage.Data);
acs800L3V = squeeze(invPowerAcs800Data.Values.L3Voltage.Data);

acs800L1U = squeeze(invPowerAcs800Data.Values.L1Current.Data);
acs800L2U = squeeze(invPowerAcs800Data.Values.L2Current.Data);
acs800L3U = squeeze(invPowerAcs800Data.Values.L3Current.Data);

acs880L1V = squeeze(invPowerAcs880Data.Values.L1Voltage.Data);
acs880L2V = squeeze(invPowerAcs880Data.Values.L2Voltage.Data);
acs880L3V = squeeze(invPowerAcs880Data.Values.L3Voltage.Data);

acs880L1U = squeeze(invPowerAcs880Data.Values.L1Current.Data);
acs880L2U = squeeze(invPowerAcs880Data.Values.L2Current.Data);
acs880L3U = squeeze(invPowerAcs880Data.Values.L3Current.Data);


shaftSpeed_rpm = squeeze(shaftSignalsData.Values.absEncoderSpeed_rpm.Data);
motorSpeed_rpm = squeeze(acs880SignalsData.Values.motorSpeed_rpm.Data);

sidAcs800TorqueSetpoint_Nm = squeeze(sidInfoSignalsData.Values.acs800TorqueSetpoint_Nm.Data);
sidAcs880SpeedSetpoint_rpm = squeeze(sidInfoSignalsData.Values.acs880SpeedSetpoint_rpm.Data);
sidFromFileCaseCounter = squeeze(sidInfoSignalsData.Values.fromFileCaseCounter.Data);

figure
subplot(2,1,1)
plot(time, expTime);
xlabel('time (s)')
ylabel('Exp. time (s)')

hold on
subplot(2,1,2)
plot(time, runCounter);
xlabel('time (s)')
ylabel('Run counter')

figure
subplot(3,1,1)
plot(time, sidFromFileCaseCounter)
xlabel('time (s)')
ylabel('SID Case conuter')

subplot(3,1,2)
plot(time, sidAcs880SpeedSetpoint_rpm)
xlabel('time (s)')
ylabel('ACS880 Speed Set (rpm)')

subplot(3,1,3)
plot(time, sidAcs800TorqueSetpoint_Nm)
xlabel('time (s)')
ylabel('ACS880 Torque Set (Nm)')



figure
subplot(2,1,1)
plot(time, acs800L1V)
hold on
plot(time, acs800L2V)
plot(time, acs800L3V)
xlabel('time (s)')
ylabel('Voltage (V)')
legend('L1','L2','L3')

subplot(2,1,2)
hold on
plot(time, acs800L1U)
hold on
plot(time, acs800L2U)
plot(time, acs800L3U)
xlabel('time (s)')
ylabel('Current (A)')
legend('L1','L2','L3')


figure
subplot(2,1,1)
plot(time, acs880L1V)
hold on
plot(time, acs880L2V)
plot(time, acs880L3V)
xlabel('time (s)')
ylabel('Voltage (V)')
legend('L1','L2','L3')

subplot(2,1,2)
hold on
plot(time, acs880L1U)
hold on
plot(time, acs880L2U)
plot(time, acs880L3U)
xlabel('time (s)')
ylabel('Current (A)')
legend('L1','L2','L3')


figure
plot(time, shaftSpeed_rpm)
hold on
plot(time, motorSpeed_rpm)
ylim([-1700 1700])
xlabel('time (s)')
ylabel('Speed (rpm)')
legend('Shaft','Motor')

return
time = hptoCtrlData.Values.runCounter.Time;
runCounter = squeeze(hptoCtrlData.Values.runCounter.Data);

excForce_N = squeeze(hptoCtrlData.Values.excForce_N.Data);

ctrlSignal1 = squeeze(hptoSignalsData.Values.ctrlSignal1.Data);
ctrlSignal2 = squeeze(hptoSignalsData.Values.ctrlSignal2.Data);
genTorqueCmd_Nm = squeeze(hptoSignalsData.Values.genTorqueCmd_Nm.Data);

% =========================================================================  

figure
subplot(4,1,1)
plot(time,runCounter)
xlabel('t (s)')
ylabel('Counter')
title('Run counter')

subplot(4,1,2)
plot(time,excForce_N/1000)
xlabel('t (s)')
ylabel('Force (kN)')
title('Excitation force')

subplot(4,1,3)
plot(time,ctrlSignal1)
hold on
plot(time,ctrlSignal2)
xlabel('t (s)')
ylabel('Ctrl signal')
title('Swash plate control')
legend('Swash plate 1','Swash plate 2')

subplot(4,1,4)
plot(time, genTorqueCmd_Nm)
xlabel('t (s)')
ylabel('Torque (Nm)')
title('Gen torque control')



