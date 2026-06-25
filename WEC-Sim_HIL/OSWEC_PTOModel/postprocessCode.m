simOutput = logsout; % Simulation output

TimeV = logsout.getElement('flowRateCV1').Values.Time;
flowRateCV1 = squeeze(logsout.getElement('flowRateCV1').Values.Data);% [lpm]
flowRateCV2 = squeeze(logsout.getElement('flowRateCV2').Values.Data);% [lpm]
flowRateCV3 = squeeze(logsout.getElement('flowRateCV3').Values.Data);% [lpm]
flowRateCV4 = squeeze(logsout.getElement('flowRateCV4').Values.Data);% [lpm]
pressureA = squeeze(logsout.getElement('pressureA').Values.Data);% [psi]
pressureB = squeeze(logsout.getElement('pressureB').Values.Data);% [psi]
pressureC = squeeze(logsout.getElement('pressureC').Values.Data);% [psi]
pressureD = squeeze(logsout.getElement('pressureD').Values.Data);% [psi]
pistonForce = squeeze(logsout.getElement('pistonForce').Values.Data);% [N]
%vel = squeeze(logsout.getElement('vel').Values.Data);% [m/s]
genTorque = squeeze(logsout.getElement('genTorque').Values.Data);% [Nm]
shaftSpeed = squeeze(logsout.getElement('shaftSpeed').Values.Data);% [rpm]
shaftPower = squeeze(logsout.getElement('shaftPower').Values.Data);% [W]
pistonPowerMech = logsout.getElement('pistonPowerMech').Values.Data;
flowRateAccHP = logsout.getElement('flowRateAccHP').Values.Data;
ptoPowerMech = logsout.getElement('ptoPowerMech').Values.Data;
%pistonVel = logsout.getElement('vel').Values.Data;
%plot(TimeV,pistonForce.*vel)


figure
subplot(3,1,1)
plot(TimeV,flowRateCV1,TimeV,flowRateCV2)
legend('Flow CV1', 'Flow CV2')

subplot(3,1,2)
plot(TimeV,pressureB,TimeV,pressureD)
legend('Pressure B', 'Pressure D')

subplot(3,1,3)
plot(TimeV,(pressureD-pressureB))
legend('DeltaP CV2')

figure
plot(TimeV,pistonPowerMech,TimeV,shaftPower)
%xlim([0 200])
legend('Input Power', 'Gen. Power')

figure
plot(TimeV,pressureC)
%xlim([0 200])
legend('Pressure HM')