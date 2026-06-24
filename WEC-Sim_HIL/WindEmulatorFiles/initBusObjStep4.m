% bus definitions for Step 4 
clearvars; close all; clc;

busFile = 'windEmulatorBusDefs.mat';

%% === ACS880 signals (emulates the hydraulic PTO incl. excitation) =======
clear acs880SignalStruct;
acs880SignalStruct.motorVoltage_V = 0.0;
acs880SignalStruct.motorCurrent_A = 0.0;
acs880SignalStruct.frequency_Hz = 0.0;
acs880SignalStruct.motorSpeed_rpm = 0.0;
acs880SignalStruct.motorTorque_Nm = 0.0;
acs880SignalStruct.shaftPower_W = 0.0;
acs880SignalStruct.statusWord = uint16(0);

acs880SignalBusInfo = Simulink.Bus.createObject(acs880SignalStruct);
acs880SignalBus = eval(acs880SignalBusInfo.busName); % assign to our prefered naming
% =========================================================================

%% === ACS880 ctrl signals ================================================
clear acs880CtrlStruct;
acs880CtrlStruct.ctrlWord = uint16(0);
acs880CtrlStruct.state = int32(0);
acs880CtrlStruct.torqueSetpoint_Nm = 0.0;
acs880CtrlStruct.torqueSetpoint_percent = 0.0;

acs880CtrlBusInfo = Simulink.Bus.createObject(acs880CtrlStruct);
acs880CtrlBus = eval(acs880CtrlBusInfo.busName); % assign to our prefered naming
% =========================================================================


%% === ACS800 signals (generator side with active front end) ===============
clear acs800SignalStruct;
acs800SignalStruct.dcBusVoltage_V = 0.0;
acs800SignalStruct.frequency_Hz = 0.0;
acs800SignalStruct.temperature = 0.0;
acs800SignalStruct.motorSpeed_rpm = 0.0;
acs800SignalStruct.motorTorque_Nm = 0.0;
acs800SignalStruct.shaftPower_W = 0.0;
acs800SignalStruct.statusWord = uint16(0);

acs800SignalBusInfo = Simulink.Bus.createObject(acs800SignalStruct);
acs800SignalBus = eval(acs800SignalBusInfo.busName); % assign to our prefered naming
% =========================================================================

%% === Shaft Signals ================================================
clear shaftSignalStruct;
shaftSignalStruct.torqueActual_Nm = 0;
shaftSignalStruct.absEncoderCounts = uint32(0);
shaftSignalStruct.absEncoderTurns = int32(0);
shaftSignalStruct.absEncoderPosition_rad = 0;
shaftSignalStruct.absEncoderSpeed_rpm = 0;
shaftSignalStruct.absEncoderStatus1 = uint8(0); % encoder has 2 uint8 status words
shaftSignalStruct.absEncoderStatus2 = uint8(0);

shaftSignalBusInfo = Simulink.Bus.createObject(shaftSignalStruct);
shaftSignalBus = eval(shaftSignalBusInfo.busName); % assign to our prefered naming
% =========================================================================

%% === ACS800 ctrl signals (could potentially be the same type as ACS880, but keep seperate for now)
clear acs800CtrlStruct;
acs800CtrlStruct.ctrlWord = uint16(0);
acs800CtrlStruct.state = int32(0);
acs800CtrlStruct.torqueSetpoint_Nm = 0.0;
acs800CtrlStruct.torqueSetpoint_percent = 0.0;

acs800CtrlBusInfo = Simulink.Bus.createObject(acs800CtrlStruct);
acs800CtrlBus = eval(acs800CtrlBusInfo.busName); % assign to our prefered naming
% =========================================================================

%% === HPTO signal bus ===================================================
clear hptoSignalStruct;
hptoSignalStruct.genTorqueCmd_Nm = 0.0;             % the setpoint for the generator torque
hptoSignalStruct.pressure_bar = 0.0;
hptoSignalStruct.hmOutputShafTorque_Nm = 0.0;       % the torque at the output (generator side) hydraulic machine
hptoSignalStruct.genShaftSpeed_rpm = 0.0;           % speed on the excitation side shaft
hptoSignalStruct.excShaftSpeed_rpm = 0.0;           % speed on the generator side shaft
hptoSignalStruct.excShaftTorque_Nm = 0.0;          % torque on the excitation side shaft
hptoSignalStruct.ctrlSignal1 = 0.0;
hptoSignalStruct.ctrlSignal2 = 0.0;
hptoSignalStruct.genPumpFlow_lpm = 0.0;

hptoSignalBusInfo = Simulink.Bus.createObject(hptoSignalStruct);
hptoSignalBus = eval(hptoSignalBusInfo.busName); % assign to our prefered naming
% =========================================================================

%% === HPTO signal bus ===================================================
clear hptoCtrlStruct;
hptoCtrlStruct.speedRef_rpm = 0.0;
hptoCtrlStruct.excForce_N = 0.0;
hptoCtrlStruct.genSpeedActual = 0.0;
hptoCtrlStruct.speedCtrlReset = boolean(false);

hptoCtrlBusInfo = Simulink.Bus.createObject(hptoCtrlStruct);
hptoCtrlBus = eval(hptoCtrlBusInfo.busName); % assign to our prefered naming
% =========================================================================


%% === SID (system characterization) ctrl signals =========================
clear sidCtrlStruct;
sidCtrlStruct.acs800Torque_Nm = 0.0;
sidCtrlStruct.acs880Torque_Nm = 0.0;

sidCtrlBusInfo = Simulink.Bus.createObject(sidCtrlStruct);
sidCtrlBus = eval(sidCtrlBusInfo.busName); % assign to our prefered naming
% =========================================================================

%% === SID info signals ===================================================
clear sidInfoStruct;
sidInfoStruct.acs800TorqueSetpoint_Nm = 0.0;
sidInfoStruct.acs880SpeedSetpoint_rpm = 0.0;
sidInfoStruct.sidType = uint32(0);
sidInfoStruct.fromFileCaseCounter = 0;

sidInfoBusInfo = Simulink.Bus.createObject(sidInfoStruct);
sidInfoBus = eval(sidInfoBusInfo.busName); % assign to our prefered naming
% =========================================================================

%% === Experiment control bus =============================================
clear expCtrlStruct;
expCtrlStruct.expType = uint16(0); % will be used as enum in code
expCtrlStruct.runHil = boolean(false);
expCtrlStruct.runSid = boolean(false);
expCtrlStruct.resetHilIntegrator = boolean(false);
expCtrlStruct.resetSidIntegrator = boolean(false);
expCtrlStruct.time = 0;
expCtrlStruct.ramp = 0;
expCtrlStruct.runCounter = uint32(0);
expCtrlStruct.stepCounter = uint32(0);

expCtrlBusInfo = Simulink.Bus.createObject(expCtrlStruct);
expCtrlBus = eval(expCtrlBusInfo.busName); % assign to our prefered naming
% =========================================================================


%% === power observations (Beckhoff) ======================================
clear invPowerStruct
invPowerStruct.L1InaccurateU = logical(false);
invPowerStruct.L1InaccurateI = logical(false);
invPowerStruct.L1Voltage = 0.0;
invPowerStruct.L1Current = 0.0;
invPowerStruct.L1PowFactor = 0.0;
invPowerStruct.L1ActivePow = 0.0;
invPowerStruct.L1THDu = 0.0;
invPowerStruct.L1THDi = 0.0;

invPowerStruct.L2InaccurateU = logical(false);
invPowerStruct.L2InaccurateI = logical(false);
invPowerStruct.L2Voltage = 0.0;
invPowerStruct.L2Current = 0.0;
invPowerStruct.L2PowFactor = 0.0;
invPowerStruct.L2ActivePow = 0.0;
invPowerStruct.L2THDu = 0.0;
invPowerStruct.L2THDi = 0.0;

invPowerStruct.L3InaccurateU = logical(false);
invPowerStruct.L3InaccurateI = logical(false);
invPowerStruct.L3Voltage = 0.0;
invPowerStruct.L3Current = 0.0;
invPowerStruct.L3PowFactor = 0.0;
invPowerStruct.L3ActivePow = 0.0;
invPowerStruct.L3THDu = 0.0;
invPowerStruct.L3THDi = 0.0;

invPowerStruct.Frequency = 0.0;
invPowerStruct.TotalPowFactor = 0.0;
invPowerStruct.TotalActivePow = 0.0;
invPowerStruct.L1L2Voltage = 0.0;
invPowerStruct.L2L3Voltage = 0.0;
invPowerStruct.L3L1Volage = 0.0;

invPowerBusInfo = Simulink.Bus.createObject(invPowerStruct);
invPowerBus = evalin('base',invPowerBusInfo.busName); % assign to our prefered naming
% =========================================================================


%% === power observations (HBM) ===========================================
clear hmbPowerStruct
hbmPowerStruct.acqTime = 0.0;               % time (s)
hbmPowerStruct.acqState = 0.0;              % state (1 for running)
hbmPowerStruct.latency = 0.0;               % latency from acquire -> send

hbmPowerStruct.driveInIRMS = 0.0;              % drive in RMS current, DC
hbmPowerStruct.driveInURMS = 0.0;              % drive in RMS voltage, DC
hbmPowerStruct.driveInIMean = 0.0;          % drive in mean current, DC
hbmPowerStruct.driveInUMean = 0.0;          % drive in mean voltage, DC
hbmPowerStruct.driveInP = 0.0;              % drive in power, DC

% powerStructHbm.driveOutCycleCheck = 0.0;    % drive out cycle check (likely the same as frequency)
% powerStructHbm.driveOutCycleCount = 0.0;    % drive out cycle count 
% powerStructHbm.driveOutCycleTimeout = 0.0;  % drive out cycle timeout
% powerStructHbm.driveOutFrequency = 0.0;     % drive out AC frequency

hbmPowerStruct.driveOutI1RMS = 0.0;          % drive out RMS current phase 1
hbmPowerStruct.driveOutI2RMS = 0.0;          % drive out RMS current phase 2
hbmPowerStruct.driveOutI3RMS = 0.0;         % drive out RMS current phase 3

hbmPowerStruct.driveOutU1RMS = 0.0;            % drive out RMS voltage phase 1
hbmPowerStruct.driveOutU2RMS = 0.0;            % drive out RMS voltage phase 2
hbmPowerStruct.driveOutU3RMS = 0.0;            % drive out RMS voltage phase 3

hbmPowerStruct.driveOutI1Mean = 0.0;            % drive out Mean current phase 1
hbmPowerStruct.driveOutI2Mean = 0.0;            % drive out Mean current phase 2
hbmPowerStruct.driveOutI3Mean = 0.0;            % drive out Mean current phase 3

hbmPowerStruct.driveOutU1Mean = 0.0;            % drive out Mean voltage phase 1
hbmPowerStruct.driveOutU2Mean = 0.0;            % drive out Mean voltage phase 2
hbmPowerStruct.driveOutU3Mean = 0.0;            % drive out Mean voltage phase 3

hbmPowerStruct.driveOutP1 = 0.0;            % drive out mean cycle power phase 1
hbmPowerStruct.driveOutP2 = 0.0;            % drive out mean cycle power phase 2
hbmPowerStruct.driveOutP3 = 0.0;            % drive out mean cycle power phase 3

% powerStructHbm.driveOutQ1 = 0.0;            % drive out reactive power phase 1 ( Q = sqrt(S^2 - P^2) )
% powerStructHbm.driveOutQ2 = 0.0;            % drive out reactive power phase 2
% powerStructHbm.driveOutQ3 = 0.0;            % drive out reactive power phase 3
% 
% powerStructHbm.driveOutS1 = 0.0;            % drive out apparent power phase 1 (RMS Current * RMS Voltage)
% powerStructHbm.driveOutS2 = 0.0;            % drive out apparent power phase 2
% powerStructHbm.driveOutS3 = 0.0;            % drive out apparent power phase 3

% powerStructHbm.driveOutIRMS = 0.0;             % drive out total RMS current (I1+I2+I3)/3  ?TODO - double check what this represents
% powerStructHbm.driveOutURMS = 0.0;             % drive out total RMS voltage (U1+U2+U3)/3  ?TODO - double check what this represents
hbmPowerStruct.driveOutP = 0.0;             % drive out total Power (P1+P2+P3)          ?TODO - double check what this represents
% powerStructHbm.driveOutQ = 0.0;             % drive out total reactive power (Q1+Q2+Q3) ?TODO - double check what this represents
% powerStructHbm.driveOutS = 0.0;             % drive out total apparent power (S1+S2+S3) ?TODO - double check what this represents

hbmPowerStruct.essOutIRMS = 0.0;               % ESS out RMS current 
hbmPowerStruct.essOutURMS = 0.0;               % ESS out RMS voltage
hbmPowerStruct.essOutIMean = 0.0;           % ESS out mean current 
hbmPowerStruct.essOutUMean = 0.0;           % ESS out mean voltage
hbmPowerStruct.essOutP = 0.0;               % ESS out cycle mean power

hbmPowerStruct.loadIRMS = 0.0;                 % Load RMS current
hbmPowerStruct.loadURMS = 0.0;                 % Load RMS voltage
hbmPowerStruct.loadIMean = 0.0;                 % Load mean current
hbmPowerStruct.loadUMean = 0.0;                 % Load mean voltage
hbmPowerStruct.loadP = 0.0;                   % Load cycle mean power

% powerStructHbm.eff1PIn = 0.0;               % Input power for efficiency calculations (same as drive in P) ?TODO - double check what these represents
% powerStructHbm.eff1POut = 0.0;              % Output power for efficiency calculations (same as drive out P) ???TODO - double check
% powerStructHbm.eff1PLossIO = 0.0;           % Power loss, with direction from input to output
% powerStructHbm.eff1PLossOI = 0.0;           % Power loss, with direction from output to input
% powerStructHbm.eff1EtaIO = 0.0;             % Efficiency from input to output                
% powerStructHbm.eff1EtaOI = 0.0;             % Efficiency from output to input

hbmPowerBusInfo = Simulink.Bus.createObject(hbmPowerStruct);
hbmPowerBus = evalin('base',hbmPowerBusInfo.busName); % assign to our prefered naming
% =========================================================================

save(busFile,   'acs800SignalStruct',...
                'acs800SignalBus',...
                'acs800CtrlStruct',...
                'acs800CtrlBus',...
                'acs880SignalStruct',...
                'acs880SignalBus',...
                'acs880CtrlStruct',...
                'acs880CtrlBus',...
                'hptoSignalStruct',...
                'hptoSignalBus',...
                'hptoCtrlStruct',...
                'hptoCtrlBus',...
                'sidCtrlStruct',...
                'sidCtrlBus',...
                'expCtrlStruct',...
                'expCtrlBus',...
                'shaftSignalStruct',...
                'shaftSignalBus',...
                'invPowerStruct',...
                'invPowerBus',...
                'hbmPowerStruct',...
                'hbmPowerBus',...
                'sidInfoStruct',...
                'sidInfoBus');
            
            
            