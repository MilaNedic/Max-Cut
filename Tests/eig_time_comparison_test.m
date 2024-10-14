% Test script for calculating the average time for computing eigenvalues
% of a random matrix of size n x n over m iterations
n = 100;
m = 10;
[cpu_time, gpu_time] = eig_time_comparison(n, m);
ratio = cpu_time/gpu_time;
fprintf('Average CPU time: %f seconds\n', cpu_time);
fprintf('Average GPU time: %f seconds\n', gpu_time);
fprintf('cpu_time/gpu_time: %f', ratio);
