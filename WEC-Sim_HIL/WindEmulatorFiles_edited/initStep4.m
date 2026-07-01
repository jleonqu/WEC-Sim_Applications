ls% Step 4 is moving to a more complete model incl. file writing and UI

clearvars; close all; clc;
addpath('icons')

% Comment from remote connection

%% === general model parameters and settings ==============================
includeHBM = boolean(false); % include or exclude the HBM power analyzer

hbmCtScale = 1/1.5; % this scale accounts for CT difference between SWEPT and DETL (1:1500 vs. 1:1000). Affects current and power.
warning('HBM scale set above. Check this when changing CTs at DETL')

mdlName = 'windEmulatorStep4';
busDefs = 'windEmulatorBusDefs.mat';
%buildDir = fullfile('c:','simulink_build');
buildDir = pwd;

currentDir = pwd;
if includeHBM
    eCATFile = 'WindTurbineEmulator-HBM-2022-12-08_4ms.xml';   % 2ms; with HBM
else
    %eCATFile = 'WindTurbineEmulator-2022-12-06_4ms.xml';    % 4ms; no HBM 
    eCATFile = 'newDevice3File_windEmulator_noHBM-2026-06-26.xml';    % 4ms; no HBM 
end

eCATFilePath = fullfile(currentDir,'ethercat',eCATFile);

appName = 'windEmulatorStep4App.mlapp';
%tgName = 'EGIBaseline2';
%tgName = 'BlueSpeedgoat';
tgName = 'speedgoat_WECHIL';

Ts = 1/250; % SLRT sample time
movingAverageLength = 10/Ts;    % number of samples for moving average (e.g. 10s)

appDecimation = 50;            % decimation for UI rate
appDecimationReduced = 250;     % reduced decimation (e.g. state indicators)
fileDecimation = 1;             % e.g. at 500Hz sample rate, 2 corresponds to 250Hz file writing

bar2pa = 1e5;
bar2psi = 14.503773773;
psi2bar = 1/bar2psi;
psi2pa  = psi2bar * 1e5;

rpm2radps = 2*pi/60;
radps2rpm = 1/rpm2radps;

m3ps2lpm = 60*1000;
lpm2m3ps = 1/m3ps2lpm;  

hp2W = 745.7; % hp to Watt

cc2m3 = (1/100)^3;

rev2rad = 2*pi; % revolutions to radians    
% =========================================================================

%% === Disable Simulink Data Inspector logging (slows down app) ===========
Simulink.sdi.setRecordData(false);
% =========================================================================

%% === initial wave conditions (can be changed in UI) =====================
% Define the wave condition
A = 1.0; % TODO JS change this to initial conditions for the app
T = 4;

waveAmpInit = A;
wavePeriodInit = T;

rampTime = 5;   % ramp time in (s) for setpoint generation

runTimeInit = 3600;     % the initial (default) run time; can be changed in UI

excAmpInit = 1000;      % initial value for force excitation amplitude in app (N)
excFreqInit = 0.25;     % initial value for force excitation frequency in app (Hz)
% =========================================================================

%% === pto specific quantities ============================================
% hydraulic motor variables
Dm_max = 40; % [cc/rev] Maximum motor displacement

% efficiency model parameters
nShaftSpeed = 2600;                         % [rpm] Nominal shaft speed
nPreDrop = 2900;                            % [psi] Nominal pressure drop
nPreDropMin = 725;                          % [psi]
EffVolNom = 0.96;                           % [1] Volumetric efficiency at nominal conditions
EffMechMax  = 0.93;                         % [1] Mechanical efficiency at nPreDrop
EffMechMin = 0.79;                          % [1] Mechanical efficiency at nPreDropMin
nTorqueNoLoad = 0.5;                        % [Nm] No torque load
FactorEff = psi2pa*cc2m3*(1/rev2rad);       % factor used to take care of the units (psi->Pa and cc/rev -> m3/rad)

% accumulator
AccVol = 1e-3;      %[m^3] Accumulator volume
Acc_PC = 1;         %[psi] Accumulator pre charge
Acc_MaxP = 5000;    %[psi] Accumulator maximum pressure

% relief pressure valve
PressureValveMax = 5000;    %[psi]
ValvePreRange = 50;         %[psi]

% generator shaft inertia
genShaftInertia_kgm2 = 0.25*1;   %[kg/m2]
% =========================================================================


%% === HPTO controller specific quantities ================================
%Torque input control
TorqueInputControl.PG = 0.05; %Proportional gain
TorqueInputControl.IG = 0.1*0; %Integral gain

%Pressure control
PressureControl.PG = 0.1*2; %Proportional gain
PressureControl.IG = 0.1*0; %Integral gain

%Speed control
SpeedControl.PG = 0.7; %Proportional gain
SpeedControl.IG = 0.3; %Integral gain - no integral gain at this point; the code holds the I term to zero via the reset currently; change this if introducing I gain here.

%Dead band speed controller
deadBandController = boolean(false);
%deadBandController = boolean(true);
deadbandTorqueSlewRate = 300;           %slew rate limit on the deadband controller
belowMinPGain = 1/200;                  %1/X means that we get full torque at X RPM below min shaft speed
ShaftSpeedRef = 1100;                   %[rpm] Max. Shaft speed reference
ShaftSpeed_min = 800;                     %[rpm] Min. shaft speed reference
TorqueLoadMax = 75;                     %[Nm] Max. torque for deadband controller

% min pressure reference
minPressureRef_psi = 1000; %[psi] Min. pressure reference % TODO - revise this limit; 1000psi leads to 43.89Nm, which is ~1/2 of the wind emulator torque

% min torque reference: Pressure(Pa) * Displacement (m3/s)  
minTorqueRef_Nm = (minPressureRef_psi*psi2pa)*(Dm_max*cc2m3*(1/rev2rad));

fprintf('A %4.2f psi pressure reference requires a minimum torque reference of %4.2f Nm\n.',minPressureRef_psi,minTorqueRef_Nm) 

% =========================================================================


%% === wavebot related quantities =========================================
load('waveBot\WaveBot2XBEM.mat'); %Load admittance model and excitation force function

% Admittance model and excitation force

Yi = WaveBot2XBEM.sysA;
Hex = WaveBot2XBEM.HexBEM;
wex = WaveBot2XBEM.wFreq;
Tex = 1./(WaveBot2XBEM.wFreq / (2*pi));
Fe = waveAmpInit*interp1(wex , Hex, (1/wavePeriodInit)*2*pi); %TODO JS - this should be done in the app, whenever we change the freq / amp

% figure
% plot(Tex, Hex/1000)
% xlabel('T (s)')
% ylabel('Fex (kN/m)')
% xlim([1 20])

% quantities to scale the gear ratio from translational to rotational 
targetOmega = 2000/60*2*pi; % target rpm on the input shaft
targetWaveBotVel = 1;     % target wavebot vel (m/s) % TODO - the waveBot excitation force for ~1m wave amplitude is quite large for our 80Nm system

wheelRadius = targetWaveBotVel/targetOmega; % for wheel and axle

fprintf('Initial wave and force conditions: \n');
fprintf('Wave amplitude:   %4.2f m\n',waveAmpInit);
fprintf('Wave period:      %4.2f s\n',wavePeriodInit);
fprintf('Excitation force: %4.2f kN\n\n',Fe/1000);
% =========================================================================

%% === damping settings ===================================================
% Load optimal values of Kpos and Kdamping:

% Use function OptimizationResultsLoad(T,A) to load the values of Kpos and
% Kdamping that correspond to the wave condition defined
% Unzip folder KposKdampingOpt.zip
OptiData = OptimizationResultsLoad(wavePeriodInit,waveAmpInit); %TODO JS - this should be done in the app, when we change T and A

kDampingInit = OptiData.Kdamping;
kSpringInit  = OptiData.Kpos;

fprintf('Initial absorption controller settings:\n');
fprintf('Kdamping = %4.2f\n', kDampingInit)
fprintf('Kspring     = %4.2f\n', kSpringInit)
%% ========================================================================


%% === ABB specific parameters ============================================

% control modes
ctrlModeSpeed = logical(false); % speed control (for control word)
ctrlModeTorque = logical(true); % torque control (for control word)

acs880TorqueSetpointScaling = 100; % this is in % (not Nm); 100% should be ~81Nm. Need to verify
acs880RatedTorque = 80.494; % Nm - TODO - check with Jon
acs880TorqueFieldbusScale = 10000;

genMaxTorque = 0.975*acs880RatedTorque;

acs880MotorVoltsScaling = 1;                                            % No scale; 1=1V
acs880MotorCurrentScaling = 100/10000;                                  % Param 46.05 (100A) / FB value (10,000)
acs880FreqScaling = 60/20000;                                           % Param 46.02 (60Hz)/ FB value (20,000)
acs880SpeedScaling = 2000/20000;                                        % Param 46.01 (2000 rpm) / FB value (20,000)
acs880TorqueScaling = acs880RatedTorque/acs880TorqueFieldbusScale;      % Rated torque / FB value (10,000) 
acs880PowerScaling = 100*hp2W/10000;                                    % Param 46.04 (100 HP) / FB value (10,000) converted to W

acs880SpeedPGainInit = 0.7;                      % p Gain for the speed PI controller
acs880SpeedIGainInit = 0.3;                      % i Gain for the speed PI controller
acs880SpeedPILimUp = 0.975*acs880RatedTorque;      % upper limit for the PI controller (Nm)
acs880SpeedPILimLo = -0.975*acs880RatedTorque;     % lower limit for the PI controller (Nm)

acs880SetpointLimUpper = 0.975*acs880RatedTorque;    % upper limit, set just before demand output (Nm)
acs880SetpointLimLower = -0.975*acs880RatedTorque;   % lower limit, set just before demand output (Nm)

acs880RateLimRising = 10000;  % slew rate of torque setpoint (Nm/s)
acs880RateLimFalling = -10000; 

%ACS800
Gen_PowerScaling = 14.9109; %kW, (1771.9 rpm)*(2*pi/60)*(80.35933 Nm) = 14.9109 kW
%Gen_REF2_Max = 100; % percent
acs800TorqueNomEng = 100;                           % (%) "engineering" unit, used for torque setpoint scaling if in torque ctrl (REF2)
acs800TorqueNomFb = 10000;                          % fieldbus equivalent for REF2
acs800SpeedNomEng = 1500;                           % (rpm) engineering unit, used for speed setpoint scaling if in speed ctrl (REF1)
acs800SpeedNomFb = 20000;                           % fieldbus equivalent for REF1
acs800TorqueRating = 80.35933;                      % N*m, 59.27 LB-FT = 80.35933 Nm

acs800DcBusVoltsScaling = 1;                        % No scale; 1=1V
acs800FreqScaling = 1/100;                          % 100 = 1Hz
acs800TempScaling = 1/10;                           % 10 = 1% (calculated IGBT temperature)
acs800SpeedScaling = 1500/20000;                    % scaling for the actual feedback line (estimated motor speed)
acs800TorqueScaling = acs800TorqueRating/10000;     % ACS800 rated torque / FB scale (10,000)
acs800PowerScaling = Gen_PowerScaling;              % TODO - confirm if this is correct; how is this defined?

%Shaft signal scales
futekTorqueScale = -1/(2^31)*10/5*100; % the negative sign aligns the Futek torque with that reported by the ABB drives
absEncoderCountsToRad = 1/2^20*2*pi; %TODO - determine encoder scaling (counts to rad)
% TODO - SET THESE GAINS VIA THE UI / startup
%set_param([mdlName,'/ACS800SpeedCtrlPGain'],'Value',num2str(eval('ACS800PGainInit')))
%set_param([mdlName,'/ACS800SpeedCtrlIGain'],'Value',num2str(eval('ACS800IGainInit')))

%set_param([mdlName,'/ACS800SpeedCtrlLimUp'],'Value',num2str(eval('ACS800PILimUp')))
%set_param([mdlName,'/ACS800SpeedCtrlLimLo'],'Value',num2str(eval('ACS800PILimLo')))
% =========================================================================

%% === type definitions ===================================================
VFD_StatusWord; % call the status word definition file

% this is the internal (Simulink StateFlow) state
Simulink.defineIntEnumType('abbStateEnum', ...
    {'undefined', ...               %00 a non state (always an error)
    'init', ...                     %01 the starting point
    'notReadyToSwitchOn',...        %02 waiting for ready and no warnings
    'readyToSwitchOn',...           %03 ready to receive makeReady command
    'operationDisabled',...         %04 waiting for enable operation
    'operationEnabled',...          %05 fully operational
    'delayOff1'},...                %06 not sure about the purpose of this state}
    0:6, ...
    'Description', 'ABB SLRT state machine', ...
	'DefaultValue', 'undefined', ...
	'HeaderFile', 'abbState.h', ...
	'DataScope', 'Exported', ...
	'AddClassNameToEnumNames', true, ...
	'StorageType', 'int32');

% enum for the ABB status word 
Simulink.defineIntEnumType('abbStatusWordEnum', ...
    {'undefined', ...           %00 
     'someStatus',...           %XX
    },...                       %TODO - define some common status words as strings
    [0,1], ...
    'Description', 'ABB Status word', ...
	'DefaultValue', 'undefined', ...
	'HeaderFile', 'abbStatusWord.h', ...
	'DataScope', 'Exported', ...
	'AddClassNameToEnumNames', true, ...
	'StorageType', 'int16');

% enum for the ABB status word 
Simulink.defineIntEnumType('abbCtrlWordEnum', ...
    {'undefined', ...           %00 
     'someCtrl',...             %3190
    },...                       %TODO - define some common ctrl words as strings
    [0,1], ...
    'Description', 'ABB Ctrl word', ...
	'DefaultValue', 'undefined', ...
	'HeaderFile', 'abbCtrlWord.h', ...
	'DataScope', 'Exported', ...
	'AddClassNameToEnumNames', true, ...
	'StorageType', 'uint16');

% enum for experiment type
Simulink.defineIntEnumType('expTypeEnum', ...
    {'off', ...     %00
     'sid',...      %01
     'hil',...      %02
    },...                      
    0:2, ...
    'Description', 'Experiment type', ...
	'DefaultValue', 'off', ...
	'HeaderFile', 'expType.h', ...
	'DataScope', 'Exported', ...
	'AddClassNameToEnumNames', true, ...
	'StorageType', 'uint16');

% enum for experiment type
Simulink.defineIntEnumType('sidTypeEnum', ...
    {'off', ...     %00
     'manual',...   %01
     'fromFile',... %02
    },...                      
    0:2, ...
    'Description', 'SID type', ...
	'DefaultValue', 'off', ...
	'HeaderFile', 'sidType.h', ...
	'DataScope', 'Exported', ...
	'AddClassNameToEnumNames', true, ...
	'StorageType', 'uint16');

 
%% === setup and compile the code =========================================

% load relevant files
load(busDefs)
load_system(mdlName)
set_param([mdlName,'/EtherCAT Init'],'config_file',eCATFilePath)

if includeHBM
    set_param([mdlName,'/readPowerHbm'],'Commented','off')
else
    set_param([mdlName,'/readPowerHbm'],'Commented','on')
end

%% === pick up the experimental data for a matrix of torque and speed =====
sidDataFile = 'experimentDataB.mat';
sidData = load(sidDataFile);

sidDataSpeed_rpm = sidData.experimentData.wShaft_rpm;
sidDataTorque_Nm = sidData.experimentData.torqueShaft_Nm;
sidDataCaseCounter = sidData.experimentData.count;
sidTime = sidData.experimentData.time;

% create timeseries for each inport
inportSpeedSignal = timeseries(sidDataSpeed_rpm, sidTime,'Name','Speed');
inportTorqueSignal = timeseries(sidDataTorque_Nm, sidTime,'Name','Torque');
inportCounterSignal = timeseries(sidDataCaseCounter, sidTime,'Name','Counter');

% specify input
inportSpeed   = [mdlName '/inportSpeed_rpm'];
inportTorque  = [mdlName '/inportTorque_Nm'];
inportCounter = [mdlName '/inportCaseCounter'];

% observed behavior:
% on: does not go to zero at end of signal -> actually extrapolates to high
% value; slow to load
%
% off: goes to zero at end of signal; very fast to load
set_param(inportSpeed,'Interpolate','off')
set_param(inportTorque,'Interpolate','off')
set_param(inportCounter,'Interpolate','off')

% configure model to use external input (the order MUST match the inport
% order)
set_param(mdlName,'LoadExternalInput','on');

inputData = Simulink.SimulationData.Dataset;
inputData = addElement(inputData,inportSpeedSignal);
inputData = addElement(inputData,inportTorqueSignal);
inputData = addElement(inputData,inportCounterSignal);

set_param(mdlName,'ExternalInput','inputData');

fromFileSpeedSlewRate = 250;   % rpm/s (both +ve & -ve)
fromFileTorqueSlewRate = 25;   % Nm/s (both +ve & -ve)
% =========================================================================


% === compile the code ====================================================
set_param(mdlName, 'RTWVerbose', 'off');
fprintf('*** Build Simulink RT (Speedgoat) code  ...\n\n')
slbuild(mdlName)

% %% === test Speedgoat connection ==========================================
% tg = slrealtime(tgName);
% try 
%    tg.connect
% %   speedgoat.setTargetTime(now,'TargetName',tgName)
% catch ME
%    fprintf('\n*** Target %s not connected. Stopping program. Check connection.\n',tgName)
%    fprintf('\n*** Matlab error \n %s \n\n',ME.getReport)   
%    return  
% end
% 
% if tg.isConnected 
%    fprintf('\n*** Target %s is connected at IP address %s. \n\n',tg.TargetSettings.name,tg.TargetSettings.address)
% end
% % =========================================================================  
% 
% %% === run the UI app =====================================================
% run(appName)
% % =========================================================================



