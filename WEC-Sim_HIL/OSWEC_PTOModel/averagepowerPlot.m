% plot_power_averages.m
% Extracts power data from a Simulink logsout dataset, calculates a 
% sliding moving time-average, and plots the results.

%% 1. Extract Data from logsout
% Retrieve the timeseries objects
ptoPower_ts = logsout.getElement('ptoPowerMech').Values;
shaftPower_ts = logsout.getElement('shaftPower').Values;
pistonPowerMech_ts = logsout.getElement('pistonPowerMech').Values;
powerHM_ts = logsout.getElement('powerHM').Values;

% Extract raw time and data arrays
time = ptoPower_ts.Time; % Assumes both signals use the same solver time
ptoPower_data = ptoPower_ts.Data;
shaftPower_data = shaftPower_ts.Data;
pistonPowerMech_data = pistonPowerMech_ts.Data;
powerHM_data = powerHM_ts.Data;

%% 2. Calculate Moving Averages
% A cumulative average from t=0 gets skewed by startup transients. 
% A sliding moving average perfectly filters out the intra-wave 
% oscillations to show the localized mean power over time.

% Define the moving average window size in seconds 
% (Tip: Set this to match your wave period, e.g., 4 to 10 seconds)
window_size_sec = 10; 

% Calculate the centered moving average using the exact time vector.
% The window looks back and forward by window_size_sec/2 to keep the 
% average perfectly aligned in phase with your data.
ptoPower_avg = movmean(ptoPower_data, [window_size_sec/2, window_size_sec/2], 'SamplePoints', time);
shaftPower_avg = movmean(shaftPower_data, [window_size_sec/2, window_size_sec/2], 'SamplePoints', time);
pistonPowerMech_avg = movmean(pistonPowerMech_data, [window_size_sec/2, window_size_sec/2], 'SamplePoints', time);
powerHM_avg = movmean(powerHM_data, [window_size_sec/2, window_size_sec/2], 'SamplePoints', time);

% Display the final steady-state averages in the Command Window
fprintf('--- Final Moving Average Results ---\n');
fprintf('ptoPowerMech Final Value: %10.2f W\n', ptoPower_avg(end));
fprintf('shaftPower Final Value:   %10.2f W\n', shaftPower_avg(end));
fprintf('pistonPowerMech Final Value:   %10.2f W\n', pistonPowerMech_avg(end));
fprintf('powerHM Final Value:   %10.2f W\n', powerHM_avg(end));

%% 3. Plot the Data
figure('Name', 'Power Moving Averages', 'Color', 'w');

% Plot the instantaneous transient signals (fainter/thinner lines)
plot(time, ptoPower_data, 'b-', 'LineWidth', 0.5, 'DisplayName', 'ptoPowerMech (Instant)');
hold on;
plot(time, shaftPower_data, 'r-', 'LineWidth', 0.5, 'DisplayName', 'shaftPower (Instant)');
hold on;
plot(time, pistonPowerMech_data, 'k-', 'LineWidth', 0.5, 'DisplayName', 'pistonPowerMech (Instant)');

% Plot the calculated moving averages as thick dashed lines
plot(time, ptoPower_avg, 'b--', 'LineWidth', 2.5, ...
    'DisplayName', sprintf('ptoPowerMech Moving Avg (%.1fs window)', window_size_sec));
    
plot(time, shaftPower_avg, 'r--', 'LineWidth', 2.5, ...
    'DisplayName', sprintf('shaftPower Moving Avg (%.1fs window)', window_size_sec));

plot(time, pistonPowerMech_avg, 'k--', 'LineWidth', 2.5, ...
    'DisplayName', sprintf('pistonPowerMech Moving Avg (%.1fs window)', window_size_sec));

hold off;

% Formatting
grid on;
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Power (W)', 'FontWeight', 'bold'); 
title('PTO Mechanical and Shaft Power: Moving Averages', 'FontWeight', 'bold');

% Add Legends
legend('Location', 'best', 'FontSize', 11);

% Optional: Auto-scale X-axis to fit the data perfectly
xlim([time(1), time(end)]);