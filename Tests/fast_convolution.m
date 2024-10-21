% Test script for fast convolution

% Initialize a random,c omplex input data and a random filter vector
data = complex(randn(4096,100),randn(4096,100));
filter = randn(16,1);

% Display CPU time
CPUtime = timeit(@()fastConvolution(data,filter));
fprintf('CPU time: %f\n', CPUtime);

% Select a GPU device
gpu = gpuDevice;
disp(gpu.Name + " GPU selected.")

% Transfer the data to the GPU
gData = gpuArray(data);
gFilter = gpuArray(filter);

% Display GPU time
GPUtime = gputimeit(@()fastConvolution(gData,gFilter));
fprintf('GPU time: %f\n', GPUtime)

CPUtimeVectorized = timeit(@()fastConvolutionVectorized(data,filter));
fprintf('CPU time vectorized: %f\n', CPUtimeVectorized);

GPUtimeVectorized = gputimeit(@()fastConvolutionVectorized(gData,gFilter));
fprintf('GPU time vectorized: %f\n', GPUtimeVectorized);

CPUspeedup = CPUtime/CPUtimeVectorized;
GPUspeedup = GPUtime/GPUtimeVectorized;

fprintf('CPU speedup: %f\n', CPUspeedup);
fprintf('GPU speedup: %f\n', GPUspeedup);

% Plot results
bar(categorical(["CPU" "GPU"]), ...
    [CPUtime CPUtimeVectorized; GPUtime GPUtimeVectorized], ...
    "grouped")
ylabel("Execution Time (s)")
legend("Unvectorized","Vectorized")

% ------------------------ AUX FUNCTIONS -------------------------------
function y = fastConvolution(data,filter)
% Zero-pad filter to the column length of data, and transform
[rows,cols] = size(data);
filter_f = fft(filter,rows);

% Create an array of zeros of the same size and class as data
y = zeros(rows,cols,'like',data);

for idx = 1:cols
    % Transform each column of data
    data_f = fft(data(:,idx));
    % Multiply each column by filter and compute inverse transform
    y(:,idx) = ifft(filter_f.*data_f);
end

end

% ----------------------------------------------------------------------

function y = fastConvolutionVectorized(data,filter)
% Zero-pad filter to the length of data, and transform
[rows,~] = size(data);
filter_f = fft(filter,rows);

% Transform each column of the input
data_f = fft(data);

% Multiply each column by filter and compute inverse transform
y = ifft(filter_f.*data_f);
end
