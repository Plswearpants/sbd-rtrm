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
kernel_size(1,:) = [70, 70];
kernel_size(2,:) = [70, 70];
kernel_size(3,:) = [50, 50];
%kernel_size(4,:) = [50, 50];
SNR = 5; 

%% Initialize as simulated kernel from TB model 
% Randomly choose n kernel slices from simulated LDoS data
load('example_data/LDoS_sim.mat');
sliceidx = randperm(size(LDoS_sim,3), num_kernels);
A0 = cell(1,num_kernels);

% Store noiseless version for comparison
A0_noiseless = cell(1,num_kernels);

% Calculate average variance and prepare noiseless kernels
avg_var = 0;
for n = 1:num_kernels
    A0_noiseless{n} = imresize(LDoS_sim(:,:,sliceidx(n)), kernel_size(n,:));
    % Need to put each slice back onto the sphere
    A0_noiseless{n} = proj2oblique(A0_noiseless{n});
    avg_var = avg_var + var(A0_noiseless{n},0,"all");
end
avg_var = avg_var / num_kernels;

% Apply universal noise based on average variance
eta_kernel = avg_var/SNR;  % Universal noise variance
for n = 1:num_kernels
    A0{n} = A0_noiseless{n} + sqrt(eta_kernel)*randn(size(A0_noiseless{n}));
    A0{n} = proj2oblique(A0{n});  % Project back to sphere after adding noise
end

% Now A0 contains a universal noise version of the noiseless kernels matching the SNR of the observation
% A0_noiseless contains the original clean kernels if needed for comparison

%% 2. Activation map generation:
%   Each pixel has probability theta of being a kernel location
theta_cap = 2*1e-4;
theta = theta_cap/2 + theta_cap/2 * rand(1, num_kernels);  % Generate num_kernels thetas capped by theta_cap

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

Y_clean = zeros(size(X0, 1:2));
for k = 1:size(X0, 3)
    Y_clean = Y_clean + convfft2(A0_noiseless{k}, X0(:,:,k));
end
    
% Add noise to match target SNR for observation
eta = var(Y_clean, 0, "all") / SNR;
Y = Y_clean + sqrt(eta) * randn(size(Y_clean));
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Control parameter panel~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load('SBD_STM_results_20241116_124411.mat');

%% Residual Quality Analysis
% Calculate residual-based quality metric
estimated_noise_var = estimate_noise(Y, 'std');
[var_ratio, residual] = computeResidualQuality(Y, Aout, Xout, estimated_noise_var);

% Generate estimated noise signal for comparison
noise_std = sqrt(estimated_noise_var);
estimated_noise = noise_std * randn(size(Y));

% Visualize residual and noise analysis
figure('Name', 'Residual Analysis');

% Residual visualization
subplot(2,2,1)
imagesc(residual); 
colorbar; 
title(sprintf('Residual Map\nvar(noise)/var(residual) = %.3f', var_ratio));
axis image;

subplot(2,2,2)
histogram(residual(:), 50, 'Normalization', 'probability');
title('Residual Distribution');
xlabel('Value'); 
ylabel('Probability');

% Noise visualization
subplot(2,2,3)
imagesc(estimated_noise);
colorbar;
title(sprintf('Estimated Noise Map\nvar(noise) = %.3e', estimated_noise_var));
axis image;

subplot(2,2,4)
histogram(estimated_noise(:), 50, 'Normalization', 'probability');
title('Estimated Noise Distribution');
xlabel('Value');
ylabel('Probability');

set(gcf, 'Position', [100, 100, 800, 600]);

%% Observation generation with controlled noise level
%SNR = ;
[Y, Y_clean, A0, num_kernels, kernel_size] = generateObservationWithSNRControl(X0, A0_noiseless, SNR);

% Estimate noise from observation (as we would do with real data)
estimated_noise_var = estimate_noise(Y, 'std');

% Calculate residual quality metric after reconstruction is done
residual_quality = computeResidualQuality(Y, Aout, Xout, estimated_noise_var);

%% Observation normalization 
rangetype = 'dynamic';  
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
dispfun{3} = @(Y, A, X, kernel_size, kplus) showims(Y,A0{3},X0(:,:,3),A,X,kernel_size,kplus,1);

% SBD settings.
initial_iteration = 3;
maxIT= 30;

params.lambda1 = [5e-2, 5e-2, 5e-2];  % regularization parameter for Phase I
params.phase2 = false;
params.kplus = ceil(0.5 * kernel_size);
params.lambda2 = [1e-2, 1e-2, 1e-2];  % FINAL reg. param. value for Phase II
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
%[Aout, Xout, bout, extras] = SBD_test_multi(Y, kernel_size, params, dispfun, A1, initial_iteration, maxIT);
[Aout, Xout, bout, extras] = SBD_test_multi_demixing(Y, kernel_size, params, dispfun, A1, initial_iteration, maxIT);

% Save the result
% Generate a unique filename for the work space
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('SBD_STM_results_%s.mat', timestamp);

% Ensure the filename is unique
counter = 1;
while exist(filename, 'file')
    counter = counter + 1;
    filename = sprintf('SBD_demixing_STM_results_3kernels%s_%d.mat', timestamp, counter);
end

% Save the specified variables to the workspace
save(filename, 'A0', 'X0', 'Aout', 'Xout', 'bout', 'extras', 'Y', 'SNR', 'A0_noiseless');

fprintf('Results saved to: %s\n', filename);

%% Visualization of Results
visualizeResults(Y, A0, Aout, X0, Xout, bout, extras);
%visualizeResults_old(Y, A0, Aout, X0, Xout, bout, extras);

