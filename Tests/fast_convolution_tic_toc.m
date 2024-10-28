%% Test script for fast convolution using tic and toc

% Initialize random, complex input data and a random filter vector
data = complex(randn(4096,100), randn(4096,100));
filter = randn(16,1);

% Number of iterations for averaging
iterations = 100;

% Initialize accumulators for CPU and GPU times
CPUtime_sum = 0;
GPUtime_sum = 0;
CPUtimeVectorized_sum = 0;
GPUtimeVectorized_sum = 0;

% Measure and average CPU time (non-vectorized)
for i = 1:iterations
    tic;
    fastConvolution(data, filter);
    CPUtime_sum = CPUtime_sum + toc;
end
CPUtime = CPUtime_sum / iterations;
fprintf('Average CPU time: %f seconds\n', CPUtime);

% Select a GPU device
gpu = gpuDevice;
disp(gpu.Name + " GPU selected.");

% Transfer the data to the GPU
gData = gpuArray(data);
gFilter = gpuArray(filter);

% Measure and average GPU time (non-vectorized)
for i = 1:iterations
    tic;
    fastConvolution(gData, gFilter);
    wait(gpu);  % Ensure that the GPU finishes the task before measuring time
    GPUtime_sum = GPUtime_sum + toc;
end
GPUtime = GPUtime_sum / iterations;
fprintf('Average GPU time: %f seconds\n', GPUtime);

% Measure and average CPU time (vectorized)
for i = 1:iterations
    tic;
    fastConvolutionVectorized(data, filter);
    CPUtimeVectorized_sum = CPUtimeVectorized_sum + toc;
end
CPUtimeVectorized = CPUtimeVectorized_sum / iterations;
fprintf('Average CPU time (vectorized): %f seconds\n', CPUtimeVectorized);

% Measure and average GPU time (vectorized)
for i = 1:iterations
    tic;
    fastConvolutionVectorized(gData, gFilter);
    wait(gpu);  % Ensure that the GPU finishes the task before measuring time
    GPUtimeVectorized_sum = GPUtimeVectorized_sum + toc;
end
GPUtimeVectorized = GPUtimeVectorized_sum / iterations;
fprintf('Average GPU time (vectorized): %f seconds\n', GPUtimeVectorized);

% Compute speedups
CPUspeedup = CPUtime / CPUtimeVectorized;
GPUspeedup = GPUtime / GPUtimeVectorized;

fprintf('CPU speedup: %f\n', CPUspeedup);
fprintf('GPU speedup: %f\n', GPUspeedup);

% Plot results
bar(categorical(["CPU" "GPU"]), ...
    [CPUtime CPUtimeVectorized; GPUtime GPUtimeVectorized], ...
    "grouped");
ylabel("Execution Time (s)");
legend("Unvectorized", "Vectorized");

% ------------------------ AUX FUNCTIONS -------------------------------
function y = fastConvolution(data, filter)
    % Zero-pad filter to the column length of data, and transform
    [rows, cols] = size(data);
    filter_f = fft(filter, rows);

    % Create an array of zeros of the same size and class as data
    y = zeros(rows, cols, 'like', data);

    for idx = 1:cols
        % Transform each column of data
        data_f = fft(data(:, idx));
        % Multiply each column by filter and compute inverse transform
        y(:, idx) = ifft(filter_f .* data_f);
    end
end

% ----------------------------------------------------------------------

function y = fastConvolutionVectorized(data, filter)
    % Zero-pad filter to the length of data, and transform
    [rows, ~] = size(data);
    filter_f = fft(filter, rows);

    % Transform each column of the input
    data_f = fft(data);

    % Multiply each column by filter and compute inverse transform
    y = ifft(filter_f .* data_f);
end
