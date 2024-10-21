function [cpu_time_mean, cpu_time_std, cpu_time_vec_mean, cpu_time_vec_std, ...
    gpu_time_mean, gpu_time_std, gpu_time_vec_mean, gpu_time_vec_std] = vectorization_comparison_tic_toc(f,n,m)
% Function for calculating the execution time of a function f
% n is the number of points in which we want to evaluate f
% m is the number of repetitions to calculate the average and standard deviation of the times

h = 10/(n-1);

% Initialize arrays to store times
cpu_time_arr = zeros(1, m);
cpu_time_vec_arr = zeros(1, m);
gpu_time_arr = zeros(1, m);
gpu_time_vec_arr = zeros(1, m);

for iter = 1:m
    % Measure CPU time - non-vectorized
    tic;
    i = 0;
    for t = 0:h:10
        i = i + 1;
        y(i) = f(t);
    end
    cpu_time_arr(iter) = toc;

    % Measure CPU time - vectorized
    tic;
    t = 0:h:10;
    y = f(t);
    cpu_time_vec_arr(iter) = toc;

    D = gpuDevice;
    % Measure GPU time - non-vectorized
    tic;
    i = 0;
    for t = gpuArray(0:h:10)  % Transfer t to the GPU
        i = i + 1;
        y(i) = gather(f(t));    % Ensure you gather the result back from GPU
    end
    gpu_time_arr(iter) = toc;
    wait(D);

    % Measure GPU time - vectorized
    tic;
    t = gpuArray(0:h:10);  % Transfer t to the GPU
    y = gather(f(t));         % Call f on the GPU array and gather results back to CPU
    gpu_time_vec_arr(iter) = toc;
    wait(D);
end

% Calculate means and standard deviations for each case
cpu_time_mean = mean(cpu_time_arr);
cpu_time_std = std(cpu_time_arr);

cpu_time_vec_mean = mean(cpu_time_vec_arr);
cpu_time_vec_std = std(cpu_time_vec_arr);

gpu_time_mean = mean(gpu_time_arr);
gpu_time_std = std(gpu_time_arr);

gpu_time_vec_mean = mean(gpu_time_vec_arr);
gpu_time_vec_std = std(gpu_time_vec_arr);

end
