function [cpu_time, cpu_time_vec, gpu_time, gpu_time_vec] = vectorization_comparison(f)
% function for calculating the execution time of a function

tic;
% Measure CPU time - non-vectorized
i = 0;
for t = 0:.01:10
    i = i + 1;
    y(i) = f(t);
end
cpu_time = toc;

tic;
% Measure CPU time - vectorized
t = 0:.01:10;
y = f(t);
cpu_time_vec = toc;

D = gpuDevice;
% Measure GPU time - non-vectorized
tic;
i = 0;
for t = gpuArray(0:.01:10)  % Transfer t to the GPU
    i = i + 1;
    y(i) = gather(f(t));    % Ensure you gather the result back from GPU
end
gpu_time = toc;
wait(D)

% Measure GPU time - vectorized
tic;
t = gpuArray(0:.01:10);  % Transfer t to the GPU
y = gather(f(t));         % Call f on the GPU array and gather results back to CPU
gpu_time_vec = toc;
wait(D)

end