clc; clear;

% Import Manopt and initialize the SBD package
run('../init_sbd');
fprintf('\n\n');

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Create synthetic multi-kernel observation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 0. Define parameters
num_kernels = 2;
n = 1;  % number of energy layers per kernel, default 1
image_size  = [300, 300];
kernel_size = zeros(num_kernels,2);
kernel_size(1,:) = [70, 70];
kernel_size(2,:) = [70, 70];
%kernel_size(3,:) = [50, 50];
%kernel_size(4,:) = [50, 50];

rangetype = 'dynamic';  
%% Initialize as simulated kernel from TB model 
% Randomly choose n kernel slices from simulated LDoS data
load('example_data/LDoS_sim.mat');
sliceidx = randperm(size(LDoS_sim,3), num_kernels);
A0 = cell(1,num_kernels);

% Store noiseless version for comparison
A0_noiseless = cell(1,num_kernels);

SNR = 10;  % Already defined above
for n = 1:num_kernels
    A0_noiseless{n} = imresize(LDoS_sim(:,:,sliceidx(n)), kernel_size(n,:));
    % Need to put each slice back onto the sphere
    A0_noiseless{n} = proj2oblique(A0_noiseless{n});
    
    % Create noisy version with same SNR as observation
    eta_kernel = var(A0_noiseless{n},0,"all")/SNR;  % noise variance for this kernel
    A0{n} = A0_noiseless{n} + sqrt(eta_kernel)*randn(size(A0_noiseless{n}));
    A0{n} = proj2oblique(A0{n});  % Project back to sphere after adding noise
end

% Now A0 contains noisy kernels matching the SNR of the observation
% A0_noiseless contains the original clean kernels if needed for comparison

%% 2. Activation map generation:
%   Each pixel has probability theta of being a kernel location
theta_cap = 2*1e-4;
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

% Initialize kernels with window option
%window_type = {};
 window_type = {'gaussian', 2};  % Example: gaussian window with alpha=2.5
% Other window options:
% window_type = 'hann';
% window_type = 'hamming';
% window_type = 'blackman';
% window_type = {'kaiser', 5};
% window_type = '';  % no window

% Initialize kernels with the selected window type
A1 = initialize_kernels(Y, num_kernels, kernel_size, kerneltype, window_type);

%% Display initialized kernels
figure;
for n = 1:num_kernels
    subplot(1, num_kernels, n);
    imagesc(A1{n});
    title(sprintf('Initial Kernel %d', n));
    colorbar;
    axis square;
end
sgtitle('Initialized Kernels');
%% 1. Settings

% A function for showing updates as RTRM runs
figure;

dispfun = cell(1,num_kernels);
dispfun{1} = @(Y, A, X, kernel_size, kplus) showims(Y,A0{1},X0(:,:,1),A,X,kernel_size,kplus,1); % here the last entry in the showims function is the energy layer index n. 
dispfun{2} = @(Y, A, X, kernel_size, kplus) showims(Y,A0{2},X0(:,:,2),A,X,kernel_size,kplus,1);
%dispfun{3} = @(Y, A, X, kernel_size, kplus) showims(Y,A0{3},X0(:,:,3),A,X,kernel_size,kplus,1);

% SBD settings.
initial_iteration = 3;
maxIT= 50;

params.lambda1 = [5e-2, 5e-2,1e-1];  % regularization parameter for Phase I
params.phase2 = true;
params.kplus = ceil(0.5 * kernel_size);
params.lambda2 = [1e-2, 1e-2, 5e-2];  % FINAL reg. param. value for Phase II
params.nrefine = 3;
params.signflip = 0.2;
params.xpos = true;
params.getbias = true;
params.Xsolve = 'FISTA';

% this is for the test phase only:
params.X0 = X0;
params.A0 = A0;

%% Run and save 
% 2. The fun part
[Aout, Xout, bout, extras] = SBD_test_multi(Y, kernel_size, params, dispfun, A1, initial_iteration, maxIT);
%[Aout, Xout, bout, extras] = SBD_test_multi_demixing_attempt(Y, kernel_size, params, dispfun, A1, initial_iteration, maxIT);

% Save the result

% Generate a unique filename for the work space
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('SBD_STM_results_%s.mat', timestamp);

% Ensure the filename is unique
counter = 1;
while exist(filename, 'file')
    counter = counter + 1;
    filename = sprintf('SBD_STM_results_3kernels%s_%d.mat', timestamp, counter);
end

% Save the specified variables to the workspace
save(filename, 'A0', 'X0', 'Aout', 'Xout', 'bout', 'extras', 'Y');

fprintf('Results saved to: %s\n', filename);

%% Visualization of the multi-kernel results
% Visualization of the multi-kernel results
figure;
num_kernels=size(A0,2);
kernel_size = zeros(num_kernels,2);
for n = 1:num_kernels
    kernel_size(n,:) = size(A0{n});
end
% Loop through each kernel
for n = 1:num_kernels
    % Create a new figure for each kernel
    figure;
    
    % Use showims to display the results for each kernel
    showims(Y, A0{n}, X0(:,:,n), Aout{n}, Xout(:,:,n), kernel_size(n,:), [], 1);
    
    % Add a title to each figure
    sgtitle(['Kernel ' num2str(n) ' Results']);
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

figure;
evaluateActivationReconstruction(X0, Xout, kernel_size, 1);