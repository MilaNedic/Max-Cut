% Data from the table
n = [10, 100, 1000, 5000, 10000]; % Input sizes
t_CPU = [0.000026, 0.006122, 0.779923, 36.279284, 172.897618]; % CPU times
t_GPU = [0.000124, 0.011438, 0.711531, 14.948885, 47.040588];  % GPU times

% Create the plot
figure;

% Plot CPU and GPU data
semilogx(n, t_CPU, '-o', 'LineWidth', 1, 'Color', [0 0 0.8], 'MarkerSize', 5, 'MarkerFaceColor', [0 0 0.8]); % CPU (blue)
hold on;
semilogx(n, t_GPU, '-o', 'LineWidth', 1, 'Color', [0.8 0 0], 'MarkerSize', 5, 'MarkerFaceColor', [0.8 0 0]); % GPU (dark red)

set(gca, 'YScale', 'log'); % Logarithmic y-axis for detailed scale

% Add labels and title
xlabel('n (Input size)', 'FontSize', 12);
ylabel('Execution time (seconds)', 'FontSize', 12);
title('Execution Times for CPU and GPU', 'FontSize', 14);

% Customize grid and appearance
grid on;

% Add legend
legend({'t_{CPU}', 't_{GPU}'}, 'Location', 'northwest', 'FontSize', 10);

% Adjust font size
set(gca, 'FontSize', 12);

hold off;
