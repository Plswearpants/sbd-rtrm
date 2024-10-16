clc; clear;

% Import Manopt and initialize the SBD package
run('../init_sbd');
fprintf('\n\n');

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Create synthetic multi-kernel observation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 0. Define parameters
num_kernels = 3;
n = 1;  % number of energy layers per kernel, default 1
image_size  = [300, 300];
kernel_size = zeros(num_kernels,2);
kernel_size(1,:) = [50, 50];
kernel_size(2,:) = [50, 50];
kernel_size(3,:) = [40, 40];

rangetype = 'dynamic';  
%% Initialize as simulated kernel from TB model 
% Randomly choose n kernel slices from simulated LDoS data
load('example_data/LDoS_sim.mat');
sliceidx = randperm(size(LDoS_sim,3), num_kernels);
A0 = cell(1,num_kernels);
for n = 1:num_kernels
    A0{n} = imresize(LDoS_sim(:,:,sliceidx(n)), kernel_size(n,:));
    % Need to put each slice back onto the sphere
    A0{n} = proj2oblique(A0{n});
end

%% 2. Activation map generation:
%   Each pixel has probability theta of being a kernel location
theta_cap = 1e-4;
theta = theta_cap/2 + theta_cap/2 * rand(1, num_kernels);  % Generate num_kernels thetas capped by theta_cap

SNR = 10;
eta = var(A0{1},0,"all")/SNR;             % additive noise variance

% GENERATE
X0 = zeros([image_size num_kernels]);
for k = 1:num_kernels
    X0_good = false;
    while ~X0_good
        X0(:,:,k) = double(rand(image_size) <= theta(k));  % activations are on / off
        X0_good = sum(X0(:,:,k) ~= 0) > 0;
    end
end

b0 = randn;

% observation generation with noise
Y = zeros(image_size);
for k = 1:num_kernels           
    Y = Y + convfft2(A0{k}, X0(:,:,k));
end

Y = Y + b0 + sqrt(eta)*randn(image_size);

%% Observation normalization 
Y = normalizeBackgroundToZeroMean3D(Y,rangetype); 
Y = proj2oblique(Y);

figure; 
imagesc(Y);
colorbar;
axis square;

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Start SBD-STM on the synthetic data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Initialize   
kerneltype = 'selected';   % existing options: 'random' or 'selected'

% Initialize kernels
A1 = initialize_kernels(Y, num_kernels, kernel_size, kerneltype);

%% 1. Settings

% A function for showing updates as RTRM runs
figure;

dispfun = cell(1,num_kernels);
dispfun{1} = @(Y, A, X, kernel_size, kplus) showims(Y,A0{1},X0(:,:,1),A,X,kernel_size,kplus,1); % here the last entry in the showims function is the energy layer index n. 
dispfun{2} = @(Y, A, X, kernel_size, kplus) showims(Y,A0{2},X0(:,:,2),A,X,kernel_size,kplus,1);
dispfun{3} = @(Y, A, X, kernel_size, kplus) showims(Y,A0{3},X0(:,:,3),A,X,kernel_size,kplus,1);
% Create a function handle for compute_kernel_quality_factors
compute_kernel_quality = cell(1,num_kernels);
compute_kernel_quality{1} = @(input_kernel) compute_kernel_quality_factors(A0{1}, input_kernel);
compute_kernel_quality{2} = @(input_kernel) compute_kernel_quality_factors(A0{2}, input_kernel);
compute_kernel_quality{3} = @(input_kernel) compute_kernel_quality_factors(A0{3}, input_kernel);


% SBD settings.
initial_iteration = 10;
maxIT= 3;

params.lambda1 = [1e-1,1e-1,1e-1];  % regularization parameter for Phase I
params.phase2 = false;
params.kplus = ceil(0.5 * k);
params.lambda2 = [5e-2, 5e-2, 5e-2];  % FINAL reg. param. value for Phase II
params.nrefine = 3;
params.signflip = 0.2;
params.xpos = true;
params.getbias = true;
params.Xsolve = 'FISTA';
params.compute_kernel_quality = compute_kernel_quality;  % Add the function handle to params

%% Run and save 
% 2. The fun part
[Aout, Xout, extras] = SBD_test_multi(Y, kernel_size, params, dispfun, A1, initial_iteration, maxIT);

% Save the result

% Generate a unique filename for the workspace
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('SBD_STM_results_%s.mat', timestamp);

% Ensure the filename is unique
counter = 1;
while exist(filename, 'file')
    counter = counter + 1;
    filename = sprintf('SBD_STM_results_3kernels%s_%d.mat', timestamp, counter);
end

% Save the specified variables to the workspace
save(filename, 'A0', 'X0', 'Aout', 'Xout', 'extras', 'Y');

fprintf('Results saved to: %s\n', filename);

%% Visualization of the multi-kernel results
% Visualization of the multi-kernel results
figure;

% Loop through each kernel
for i = 1:num_kernels
    % Create a new figure for each kernel
    figure;
    
    % Use showims to display the results for each kernel
    showims(Y, A0{i}, X0(:,:,i), Aout{i}, Xout(:,:,i), kernel_size(i,:), [], 1);
    
    % Add a title to each figure
    sgtitle(['Kernel ' num2str(i) ' Results']);
end

% Display the original image and the reconstructed image
figure;

% Original image
subplot(121);
imagesc(Y(:,:,1));
title('Original Image');
colorbar;
axis image;

% Reconstructed image
Y_reconstructed = zeros(size(Y));
for i = 1:num_kernels
    Y_reconstructed = Y_reconstructed + convfft2(Aout{i}, Xout(:,:,i));
end

subplot(122);
imagesc(Y_reconstructed(:,:,1));
title('Reconstructed Image');
colorbar;
axis image;

% Add a main title
sgtitle('Original vs Reconstructed Image');

% Display the difference image
figure;
imagesc(Y(:,:,1) - Y_reconstructed(:,:,1));
title('Difference Image (Original - Reconstructed)');
colorbar;
axis image;

% Print out some quantitative metrics
mse = mean((Y(:,:,1) - Y_reconstructed(:,:,1)).^2, 'all');
psnr = 10 * log10(max(Y(:,:,1), [], 'all')^2 / mse);

fprintf('Mean Squared Error: %f\n', mse);
fprintf('Peak Signal-to-Noise Ratio: %f dB\n', psnr);
