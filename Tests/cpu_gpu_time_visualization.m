% Script for visualizing results
% CPU Non-Vectorized Times (in seconds)
cpu_non_vectorized_times = [0.000006, 0.000099, 0.000180, 0.000674];

% CPU Vectorized Times (in seconds)
cpu_vectorized_times = [0.000003, 0.000008, 0.000039, 0.000231];

% GPU Non-Vectorized Times (in seconds)
gpu_non_vectorized_times = [0.000760, 0.005839, 0.056800, 0.577536];

% GPU Vectorized Times (in seconds)
gpu_vectorized_times = [0.000087, 0.000082, 0.000101, 0.000202];

% Input size 
n = [11, 101, 1001, 10001];

% ------------------------------- PLOT -----------------------------------
figure;
hold on;
plot(n, cpu_non_vectorized_times, '-o', 'Color', [0 0 0.5], 'MarkerFaceColor', [0 0 0.5], 'DisplayName', 'CPU Non-Vectorized'); 
plot(n, cpu_vectorized_times, '-o', 'Color', [0.6 0.8 1], 'MarkerFaceColor', [0.6 0.8 1], 'DisplayName', 'CPU Vectorized'); 
plot(n, gpu_non_vectorized_times, '-o', 'Color', [0.5 0 0], 'MarkerFaceColor', [0.5 0 0], 'DisplayName', 'GPU Non-Vectorized'); 
plot(n, gpu_vectorized_times, '-o', 'Color', [1 0.6 0.6], 'MarkerFaceColor', [1 0.6 0.6], 'DisplayName', 'GPU Vectorized'); 
% Setting labels and title
xlabel('n (Input size)');
ylabel('Execution time (seconds)');
title('Execution Times for CPU and GPU (Vectorized vs Non-Vectorized)');

% Log scale for y-axis
set(gca, 'YScale', 'log');

legend('show');
grid on;
hold off;

