% Test script for testing the difference in execution time of non-vectorized
% and vectorized operation on CPU and GPU

% Input parameters
f = @(x) sin(x);  
n = 11;         
m = 100;           

[cpu_time_mean, cpu_time_std, cpu_time_vec_mean, cpu_time_vec_std, ...
    gpu_time_mean, gpu_time_std, gpu_time_vec_mean, gpu_time_vec_std] = vectorization_comparison_tic_toc(f, n, m);

% Display CPU non-vectorized time (mean and standard deviation)
fprintf('CPU non-vectorized time: %f seconds (mean), %f seconds (std)\n', round(cpu_time_mean, 8), round(cpu_time_std, 8));

% Display CPU vectorized time (mean and standard deviation)
fprintf('CPU vectorized time: %f seconds (mean), %f seconds (std)\n', round(cpu_time_vec_mean, 8), round(cpu_time_vec_std, 8));

% Display GPU non-vectorized time (mean and standard deviation)
fprintf('GPU non-vectorized time: %f seconds (mean), %f seconds (std)\n', round(gpu_time_mean, 8), round(gpu_time_std, 8));

% Display GPU vectorized time (mean and standard deviation)
fprintf('GPU vectorized time: %f seconds (mean), %f seconds (std)\n', round(gpu_time_vec_mean, 8), round(gpu_time_vec_std, 8));

% Calculate speedups
cpu_speedup = cpu_time_mean / cpu_time_vec_mean;
gpu_speedup = gpu_time_mean / gpu_time_vec_mean;

% Display CPU and GPU speedup
fprintf('CPU speedup: %f\n', cpu_speedup);
fprintf('GPU speedup: %f\n', gpu_speedup);

% Display CPU vectorized time / GPU vectorized time ratio
fprintf('Average CPU vectorized time / Average GPU vectorized time: %f\n', cpu_time_vec_mean / gpu_time_vec_mean);
