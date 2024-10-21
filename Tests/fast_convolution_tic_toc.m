% Test script for fast convolution using tic and toc

% Initialize random, complex input data and a random filter vector
data = complex(randn(4096,100), randn(4096,100));
filter = randn(16,1);

% Measure CPU time (non-vectorized) using tic/toc
tic;
fastConvolution(data, filter);
CPUtime = toc;
fprintf('CPU time: %f seconds\n', CPUtime);

% Select a GPU device
gpu = gpuDevice;
disp(gpu.Name + " GPU selected.")

% Transfer the data to the GPU
gData = gpuArray(data);
gFilter = gpuArray(filter);

% Measure GPU time (non-vectorized) using tic/toc
tic;
fastConvolution(gData, gFilter);
wait(gpu);  % Ensure that the GPU finishes the task before measuring time
GPUtime = toc;
fprintf('GPU time: %f seconds\n', GPUtime);

% Measure CPU time (vectorized) using tic/toc
tic;
fastConvolutionVectorized(data, filter);
CPUtimeVectorized = toc;
fprintf('CPU time (vectorized): %f seconds\n', CPUtimeVectorized);

% Measure GPU time (vectorized) using tic/toc
tic;
fastConvolutionVectorized(gData, gFilter);
wait(gpu);  % Ensure that the GPU finishes the task before measuring time
GPUtimeVectorized = toc;
fprintf('GPU time (vectorized): %f seconds\n', GPUtimeVectorized);

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
