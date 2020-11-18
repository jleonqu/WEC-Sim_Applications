%Run bemio in the hydrodata folder and then run this script

wecSim
%Constants
topA = 0.0378;
botA = 0.0378;
pistonStroke = 10;

SimResults_RM3PTOSim = cell(1);%Original results for incompressible case using PTO-Sim
SimResults_JustPiston = cell(1);%Simulation of the piston with incompressible

% load_system('RM3_Hydraulic_PTO.slx')
% SimResults_RM3PTOSim = sim('RM3_Hydraulic_PTO.slx');

load_system('PistonModelComparison.slx')
SimResults_JustPiston = sim('PistonModelComparison.slx');

figure(7)
subplot(3,1,1)
plot(PistonVelocityOriginal(:,1),PistonVelocityOriginal(:,2))
legend('Velocity Input')
ylabel('Velocity [m/s]')

subplot(3,1,2)
plot(PistonForceOriginal(:,1),PistonForceOriginal(:,2), SimResults_JustPiston(1).ForcePistonPTOSim(:,1),SimResults_JustPiston(1).ForcePistonPTOSim(:,2),SimResults_JustPiston(1).ForceASimscapePiston(:,1),SimResults_JustPiston(1).ForceASimscapePiston(:,2),SimResults_JustPiston(1).ForceBSimscapePiston(:,1),SimResults_JustPiston(1).ForceBSimscapePiston(:,2))
legend('Original Simulation','PTO-Sim Piston','Simscape Piston A','Simscape Piston B')
ylabel('Force [N]')

subplot(3,1,3)
plot(PistonVolFlowOriginal(:,1),PistonVolFlowOriginal(:,2), SimResults_JustPiston(1).volFlowPTOSim(:,1),SimResults_JustPiston(1).volFlowPTOSim(:,2),SimResults_JustPiston(1).FlowA_Simscape(:,1),SimResults_JustPiston(1).FlowA_Simscape(:,2))
legend('Original Simulation','PTO-Sim Piston','Simscape Piston')
ylabel('Flow [m^3/s]')
xlabel('Time [s]')

% subplot(4,1,3)
% plot(PistonVolFlowOriginal(:,1),PistonVolFlowOriginal(:,2), SimResults_JustPiston(1).volFlowPTOSim(:,1),SimResults_JustPiston(1).volFlowPTOSim(:,2),SimResults_JustPiston(1).FlowA_Simscape(:,1),SimResults_JustPiston(1).FlowA_Simscape(:,2))
% legend('Original Simulation','PTO-Sim Piston','Simscape Piston')
% ylabel('Flow [m^3/s]')
% xlabel('Time [s]')