% Wave parameters
%A = 2; %[m] wave amplitude
%Testing commit
T = 4; %[s] wave period

%load('WaveBotAdmittanceModel\WaveBot3XBEM.mat'); %Load admittance model and excitation force function

%Yi = WaveBot3XBEM.sysA;
%Hex = WaveBot3XBEM.HexBEM;
%SampleTime = 0.01;
f0 = 1/T;
%Fe = A*interp1(WaveBot3XBEM.wFreq, WaveBot3XBEM.HexBEM, f0*2*pi);


%% Hydraulic system parameters
%conversion factors
in2m = 0.0254;
bar2psi = 1e5/6894.75;
pa2psi = 1/6894.75;
psi2bar = 1/bar2psi;
radsec2rpm = 60/(2*pi);
m3persecond2lpm = 1000*60;

% Hydraulic Cylinder
Dp_out = 16; %[in] External piston diameter
Dp_in = 2; %[in] Piston rod diameter
areaHC = (0.25*pi*(Dp_out*in2m)^2 - 0.25*pi*(Dp_in*in2m)^2); %[m^2] Hydraulic cylinder area
strokePiston = 5; %[m] Piston stroke
deadVolume = strokePiston*areaHC*0.01; %[m^3]

%Check valves
pCrack = 10; %[psi] Cracking Pressure
pMaxValve = 25; %[psi] Pressure for max. area
areaMaxValve = 0.25*6.02e-4; %[m^2] based on orifice equation
%areaMaxValve = 8.51e-4; %[m^2] based on orifice equation

pCrackLP = 10; %[psi] Cracking Pressure
%pMaxValveLP = 75; %[psi] Pressure for max. area
pMaxValveLP = 25; %[psi] Pressure for max. area
areaMaxValveLP = 0.25*6.02e-4; %[m^2] based on orifice equation
%areaMaxValveLP = 7.5e-4; %[m^2] based on orifice equation
%areaMaxValveLP = 8.516e-4; %[m^2] based on orifice equation

%High pressure Hydraulic accumulator
accVolHP = 50; %[liter]
pPreLoadHP = 500; %[psi] preload pressure
pMaxAccHP = 5000; %[psi]

%Low pressure Hydraulic accumulator
accVolLP = 50; %[liter]
pPreLoadLP = 100; %[psi] preload pressure
pMaxAccLP = 1000; %[psi]

%Hydraulic Motor
Dmax = 15; %[cc/rev] Max. displacement
shaftSpeedNominal = 1000; %[rpm] Nom. shaft speed
pressureNom = 4000; %[psi] Nominal pressure
kinematicVisNom = 18; %[cst] Nominal Kinematic viscosity
fluidDensityNom = 900; %[kg/m^3] Nominal fluid density
volEffNom = 0.92; % Vol. eff at nominal conditions
noLoadTorque = 0.0005; %[Nm] No load Torque
frictionTorqueVsPressure = 6e-8; %[N*m/Pa] Friction torque vs. pressure coefficient

% Crank variables
crankLength = 3;
offsetLength = 1.3;
rodLength = 5;

% Electric generator
genShaftSpeedRef = 1800; % [rpm]
pGainGen = 5;
iGainGen = 1;