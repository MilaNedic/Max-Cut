% Test script for testing the difference in execution time of non-vectorized
% and vectorized operation on CPU and GPU
f = @(x) sin(x);

[cpu_time, cpu_time_vec, gpu_time, gpu_time_vec] = vectorization_comparison(f);

fprintf('CPU non-vectorized time: %f seconds\n', round(cpu_time, 8));
fprintf('CPU vectorized time: %f seconds\n', round(cpu_time_vec, 8));
fprintf('GPU non-vectorized time: %f seconds\n', round(gpu_time, 8));
fprintf('GPU vectorized time: %f seconds\n', round(gpu_time_vec, 8));

cpu_speedup = cpu_time/cpu_time_vec;
gpu_speedup = gpu_time/gpu_time_vec;

fprintf('CPU speedup: %f\n', cpu_speedup);
fprintf('GPU speedup: %f\n', gpu_speedup);
