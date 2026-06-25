% plot_steady_state_pto.m
% This script extracts angular velocity and hydraulic PTO torque from 
% a Simulink logsout dataset, crops it to the last 100 seconds, 
% normalizes the data to a [-1, 1] range, and plots them together.

%% 1. Extract Data from logsout
% Retrieve the timeseries objects from the logsout dataset
vel_ts = logsout.getElement('angVelocity').Values;
torque_ts = logsout.getElement('ptoTorqueHydraulic').Values;
pistonForce_ts = logsout.getElement('pistonForce').Values;

% Extract raw time and data arrays
% (Assuming both signals share the same time vector from the solver)
time = vel_ts.Time;
vel_data = vel_ts.Data;
torque_data = torque_ts.Data;
pistonForce_data = pistonForce_ts.Data;

%% 2. Crop to the Last 100 Seconds (Steady State)
t_end = time(end);
t_start = max(0, t_end - 100); % Ensure it doesn't try to look past t=0

% Create a logical index for the last 100 seconds
idx = time >= t_start & time <= t_end;

% Apply the index to crop the arrays
time_steady = time(idx);
vel_steady = vel_data(idx);
torque_steady = torque_data(idx);
pistonForce_steady = pistonForce_data(idx);

%% 3. Peak Normalization
% Normalize by the maximum absolute value in the cropped window.
% This preserves zero-crossings while fitting data into the [-1, 1] range.
vel_norm = vel_steady / max(abs(vel_steady));
torque_norm = torque_steady / max(abs(torque_steady));
pistonForce_norm = pistonForce_steady / max(abs(pistonForce_steady));

%% 4. Plot the Normalized Data
figure('Name', 'Steady State PTO Dynamics', 'Color', 'w');

% Plot both lines
plot(time_steady, vel_norm, 'b-', 'LineWidth', 1.5);
hold on;
plot(time_steady, torque_norm, 'r--', 'LineWidth', 1.5);

% Plot a zero-reference line to easily spot sign changes
yline(0, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off'); 
hold off;

% Formatting
grid on;
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Normalized Amplitude [-1 to 1]', 'FontWeight', 'bold');
title('Steady State: Angular Velocity vs. Hydraulic Torque', 'FontWeight', 'bold');

% Add Legends
legend('Angular Velocity', 'PTO Hydraulic Torque', ...
       'Location', 'best', 'FontSize', 11);

% Optional: Set axis limits to make it look clean
xlim([t_start, t_end]);
ylim([-1.1, 1.1]);

%% Force plot
figure('Name', 'Steady State PTO Dynamics', 'Color', 'w');

% Plot both lines
plot(time_steady, vel_norm, 'b-', 'LineWidth', 1.5);
hold on;
plot(time_steady, pistonForce_norm, 'r--', 'LineWidth', 1.5);

% Plot a zero-reference line to easily spot sign changes
yline(0, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off'); 
hold off;

% Formatting
grid on;
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Normalized Amplitude [-1 to 1]', 'FontWeight', 'bold');
title('Steady State: Angular Velocity vs. Hydraulic Piston Force', 'FontWeight', 'bold');

% Add Legends
legend('Angular Velocity', 'PTO Hydraulic Force', ...
       'Location', 'best', 'FontSize', 11);

% Optional: Set axis limits to make it look clean
xlim([t_start, t_end]);
ylim([-1.1, 1.1]);

%% Torque and force plot
figure('Name', 'Steady State PTO Dynamics', 'Color', 'w');

% Plot both lines
plot(time_steady, torque_norm, 'b-', 'LineWidth', 1.5);
hold on;
plot(time_steady, pistonForce_norm, 'r--', 'LineWidth', 1.5);

% Plot a zero-reference line to easily spot sign changes
yline(0, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off'); 
hold off;

% Formatting
grid on;
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Normalized Amplitude [-1 to 1]', 'FontWeight', 'bold');
title('Steady State: torque vs. Force', 'FontWeight', 'bold');

% Add Legends
legend('Torque', 'Force', ...
       'Location', 'best', 'FontSize', 11);

% Optional: Set axis limits to make it look clean
xlim([t_start, t_end]);
ylim([-1.1, 1.1]);

%% Force, Torque and Velocity
figure('Name', 'Steady State PTO Dynamics', 'Color', 'w');

% Plot both lines
plot(time_steady, torque_norm, 'b-', 'LineWidth', 1.5);
hold on;
plot(time_steady, pistonForce_norm, 'r-', 'LineWidth', 1.5);
hold on;
plot(time_steady, vel_norm, 'g-', 'LineWidth', 1.5);

% Plot a zero-reference line to easily spot sign changes
yline(0, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off'); 
hold off;

% Formatting
grid on;
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Normalized Amplitude [-1 to 1]', 'FontWeight', 'bold');
title('Steady State: Torque, Force and Velocity', 'FontWeight', 'bold');

% Add Legends
legend('Torque', 'Force', 'Velocity', ...
       'Location', 'best', 'FontSize', 11);

% Optional: Set axis limits to make it look clean
%xlim([380, t_end]);
xlim([380, 384]);
ylim([-1.1, 1.1]);