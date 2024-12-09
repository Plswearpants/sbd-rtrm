clc; clear;
run('../init_sbd');

%% Load pre-generated synthetic datasets
load('results\parallel results\synthetic_datasets_20241127_170302/synthetic_datasets_20241127_170302.mat');  

%% Initialize kernels for all datasets
A1_all = initialize_all_kernels(datasets);

%% Define system parameters and create configs
lambda1_range = [0.3, 0.5];
mini_loop_range = [1, 3, 9];
maxIT = 30;

% Create parameter combinations
[lambda1_grid, mini_loop_grid] = meshgrid(lambda1_range, mini_loop_range);
param_combinations = [lambda1_grid(:), mini_loop_grid(:)];
num_param_combos = size(param_combinations, 1);

% Create parameter-specific configs upfront
root_dir = fileparts(pwd);
config_dir = fullfile(root_dir, 'config', 'worker_configs');  % Base config directory
if ~exist(config_dir, 'dir')
    mkdir(config_dir);
end

% Create config for each parameter combination
for m = 1:num_param_combos
    param_dir = fullfile(config_dir, sprintf('config_param_%d', m));
    if ~exist(param_dir, 'dir')
        mkdir(param_dir);
    end
    
    % Copy and update config files
    copyfile(fullfile(root_dir, 'config', 'Xsolve_config.mat'), ...
        fullfile(param_dir, 'Xsolve_config_tunable.mat'));
    copyfile(fullfile(root_dir, 'config', 'Asolve_config.mat'), ...
        fullfile(param_dir, 'Asolve_config_tunable.mat'));
    
    % Update configs with current parameters
    update_config(fullfile(param_dir, 'Xsolve_config_tunable.mat'), ...
        'MAXIT', param_combinations(m,2));
    update_config(fullfile(param_dir, 'Asolve_config_tunable.mat'), ...
        'options.maxiter', param_combinations(m,2));
end

%% Run parallel tests
if isempty(gcp('nocreate'))
    parpool(6);
end

% Create results directory if it doesnt exist
results_dir = 'results/parallel_tests_v2';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

for n = 1:numel(datasets)   
    fprintf('Starting Dataset %d\n', n);
    
    % Create local copies of all needed variables
    local_Y = datasets(n).Y;
    local_kernel_size = datasets(n).params.kernel_size;
    local_X0 = datasets(n).X0;
    local_A0 = datasets(n).A0;
    local_A1 = A1_all{n};
    local_param_combinations = param_combinations;  % Make local copy
    local_maxIT = maxIT;  % Make local copy
    local_dataset_num = n;  % Add this line
    
    % set parameters
    local_params = struct();
    local_params.phase2 = false;
    local_params.X0 = local_X0;
    local_params.A0 = local_A0;
    local_params.Xsolve = 'FISTA';
    local_params.xpos = true;
    local_params.getbias = true;
    local_params.debug = true;
    local_params.dataset_num = local_dataset_num;  % Add this line

    parfor param_idx = 1:num_param_combos
        param_dir = fullfile(config_dir, sprintf('config_param_%d', param_idx));
        try
            % Add input parameter validation
            fprintf('\nValidating inputs for dataset %d, param combo %d:\n', n, param_idx);
            fprintf('Y size: %s\n', mat2str(size(local_Y)));
            fprintf('kernel_size: %s\n', mat2str(local_kernel_size));
            fprintf('X0 size: %s\n', mat2str(size(local_X0)));
            fprintf('A0 size: %s\n', mat2str(size(local_A0{1})));  % Assuming cell array
            fprintf('A1 size: %s\n', mat2str(size(local_A1{1})));  % Assuming cell array
            
            % Add debug prints
            fprintf('Processing param_idx %d with lambda1=%.3f, mini_loop=%d\n', ...
                param_idx, local_param_combinations(param_idx,1), local_param_combinations(param_idx,2));
            
            [Aout, Xout, bout, extras] = SBD_test_multi_parallel(...
                local_Y, ...
                local_kernel_size, ...
                local_params, ...
                local_A1, ...
                local_param_combinations,...
                param_idx, ...
                local_maxIT);

            % Generate filename with timestamp and save in results directory
            timestamp = datestr(now, 'yyyymmdd_HHMMSS');
            filename = fullfile(results_dir, ...
                sprintf('SBD_parallel_v2_dataset%d_param_idx%d_%s.mat', ...
                n, param_idx, timestamp));

            % Create a scalar structure (1x1) with all variables to save
            s = struct();  % Initialize empty scalar struct
            s(1).Aout = {Aout};  % Explicitly wrap in cell
            s(1).Xout = {Xout};
            s(1).bout = {bout};
            s(1).extras = {extras};
            s(1).param_combinations = local_param_combinations;
            s(1).param_idx = param_idx;
            s(1).dataset_A0 = {local_A0};
            s(1).dataset_X0 = {local_X0};

            % Save results immediately using -fromstruct option
            save(filename, '-fromstruct', s);
            
            fprintf('Results saved for dataset %d, param combo %d to: %s\n', ...
                n, param_idx, filename);

        catch ME
            % Enhanced error reporting
            fprintf('Error Details for dataset %d, param combo %d:\n', n, param_idx);
            fprintf('Error Message: %s\n', ME.message);
            fprintf('Error Stack:\n');
            for k = 1:length(ME.stack)
                fprintf('  File: %s\n  Line: %d\n  Function: %s\n', ...
                    ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
            end
        end
    end
end

%% Cleanup worker configs after completion
rmdir(config_dir, 's');

fprintf('All processing complete!\n'); 