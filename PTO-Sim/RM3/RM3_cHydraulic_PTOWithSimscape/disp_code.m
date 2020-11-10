delta_P_vector = linspace(0,25e6,1000);
disp_vector = zeros(1,length(delta_P_vector));

for i = 1:1:length(delta_P_vector)
    deltaP_i = delta_P_vector(i);
    disp_vector(i) = variableMotorVolume(deltaP_i)*(100^3)*2*pi;
end

figure(1)
plot(delta_P_vector, disp_vector)