function [cpu_time, gpu_time] = eig_time_comparison(n, m)
    % function for calculating the average time of computation of
    % eigenvalues for randomly generated matrices of size n x n
    % The average is taken over m iterations for both CPU and GPU
    
    % Preallocate arrays for times
    cpu_times = zeros(1, m);
    gpu_times = zeros(1, m);
    
    for i = 1:m
        % Measure CPU time
        A = rand(n); % generate a random matrix of size n x n
        tic;
        eig(A); % compute eigenvalues
        cpu_times(i) = toc;
        
        D = gpuDevice;
        % Measure GPU time
        B = gpuArray(A); % B is a random matrix of size n x n on the GPU
        wait(D); % Ensure all previous GPU operations are done before starting timer
        tic;
        eig(B); % Compute eigenvalues on GPU
        wait(D); % Wait for GPU to complete the eigenvalue computation
        gpu_times(i) = toc; %gputimeit(@() eig(B));

    end
       
    % Take the average time of m calculations and round to 8 decimal places
    cpu_time = round(mean(cpu_times), 8);
    gpu_time = round(mean(gpu_times), 8);
end
