% init script for the wind emulator

clearvars; close all; clc;

VariableDefinition

%% === model parameters and settings ======================================
mdlName = 'windEmulatorStep3';
%eCATFile = 'D:\src\sandia\ptoTest\WindEmulatorFiles\ethercat\WindTurbineEmulator-03-24-2022.xml'; % JS environment
eCATFile = 'C:\Users\SWEPT\Desktop\Jorge\WIndTurbineEmulator\WindEmulatorFiles\ethercat\WindTurbineEmulator-03-24-2022.xml'; % Jorge environment

ctrlModeSpeed = logical(false);
ctrlModeTorque = logical(true);

rampTime = 5;   % ramp time in (s) for setpoint generation

ACS800PGainInit = 0.6; % p Gain for the speed PI controller
ACS800IGainInit = 0.0;    % i Gain for the speed PI controller
ACS800PILimUp = 80;     % upper limit in % of torque
ACS800PILimLo = -80;     % lower limit in % of torque

ACS880PILimUp = 80;     % upper limit in % of torque
ACS880PILimLo = -80;     % lower limit in % of torque

% =========================================================================

%% === wavebot related quantities =========================================
WaveBotInfo = open('waveBot/WaveBot_Zi_Yi_v_enc_F_des.mat');
Yi = WaveBotInfo.Yi;

targetOmega = 2000/60*2*pi; % target 1000 rpm on the input shaft
targetWaveBotVel = 2.5;     % target 2.5m/s wavebot vel

wheelRadius = targetWaveBotVel/targetOmega; % for wheel and axle

ratedTorque = 80.494; % Nm
genMaxTorque = 0.8*ratedTorque;
% =========================================================================


%% === setup the various scaling parameters ===============================

%ACS880
Motor_SpeedScaling = 2000; %rpm, ABB ACS880 Parameter 46.01
Motor_CurrentScaling = 100; %A
Motor_FrequencyScaling = 60; %Hz
Motor_PowerScaling = 74.5699872; %kW, 100 hp = 74.5699872 kW
Motor_TorqueScaling = 100; % this is in % (not Nm); 100% should be ~81Nm. Need to verify
TorqueFieldbusScale = 10000;

%ACS800
Gen_PowerScaling = 14.9109; %kW, (1771.9 rpm)*(2*pi/60)*(80.35933 Nm) = 14.9109 kW
Gen_REF2_Max = 100; % percent
Gen_SpeedScaling = 1500; %rpm, must be maximum of abs(Par20.01) and Par20.02
speedScalingACS800 = 1500/20000; % scaling for the actual feedback line (estimated motor speed)
Gen_TorqueScaling = 80.35933; %N*m, 59.27 LB-FT = 80.35933 Nm
%Line_ScaleACT1 = 1; % currently unused
%Line_ScaleACT2 = 1; % currently unused
% =========================================================================

%% === type definitions ===================================================
VFD_StatusWord; % call the status word definition file


%% === setup and compile the code =========================================
load_system(mdlName)
set_param([mdlName,'/EtherCAT Init'],'config_file',eCATFile)

% setup the various ACS880 switches to their default position
set_param([mdlName,'/ASC880SourceSwitch'],'sw','0') % sets to const or '0', or sine wave or '1'
set_param([mdlName,'/ACS880RampEnable'],'sw','0')
set_param([mdlName,'/ASC880MakeReady'],'sw','0')
set_param([mdlName,'/ASC880EnableOperation'],'sw','0')
set_param([mdlName,'/ASC880ResetFault'],'sw','0')
set_param([mdlName,'/ASC880EnableManualCW'],'sw','0')
set_param([mdlName,'/ACS880CtrlMode'],'value','ctrlModeTorque') % use the string rather than the value to be more explicit on the diagram

% setup the various ACS800 switches to their default position
set_param([mdlName,'/ASC800SourceSwitch'],'sw','0') % sets to const or '0', or sine wave or '1'
set_param([mdlName,'/ACS800RampEnable'],'sw','0')
set_param([mdlName,'/ASC800MakeReady'],'sw','0')
set_param([mdlName,'/ASC800EnableOperation'],'sw','0')
set_param([mdlName,'/ASC800ResetFault'],'sw','0')
set_param([mdlName,'/ASC800EnableManualCW'],'sw','0')
set_param([mdlName,'/ACS800CtrlMode'],'value','ctrlModeTorque') % use the string rather than the value to be more explicit on the diagram
set_param([mdlName,'/ACS800ResetIntegrator'],'sw','0')

set_param([mdlName,'/ACS800SpeedCtrlPGain'],'Value',num2str(eval('ACS800PGainInit')))
set_param([mdlName,'/ACS800SpeedCtrlIGain'],'Value',num2str(eval('ACS800IGainInit')))

set_param([mdlName,'/ACS800SpeedCtrlLimUp'],'Value',num2str(eval('ACS800PILimUp')))
set_param([mdlName,'/ACS800SpeedCtrlLimLo'],'Value',num2str(eval('ACS800PILimLo')))

% compile the code
set_param(mdlName, 'RTWVerbose', 'off');
fprintf('*** Build Simulink RT (Speedgoat) code  ...\n\n')
slbuild(mdlName)
open_system(mdlName)
