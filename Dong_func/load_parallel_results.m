% Load and organize parallel test results for visualization
function [metrics] = load_parallel_results()
    % Initialize storage structure
    metrics = struct();
    
    % Get folder path from user
    default_path = fullfile(pwd, 'examples');
    folder_path = uigetdir(default_path, 'Select folder containing SBD parallel results');
    
    if folder_path == 0  % User canceled
        error('Folder selection canceled by user');
    end
    
    % Get list of all result files
    files = dir(fullfile(folder_path, 'SBD_parallel_dataset*.mat'));
    
    if isempty(files)
        error('No SBD parallel result files found in selected folder');
    end
    
    % Display loading information
    fprintf('Found %d result files in:\n%s\n\n', length(files), folder_path);
    
    % Extract unique dataset numbers and parameter indices
    dataset_nums = [];
    param_indices = [];
    for i = 1:length(files)
        parts = split(files(i).name, {'_', '.'});
        dataset_nums = [dataset_nums; str2double(parts{3}(8:end))];
        param_indices = [param_indices; str2double(parts{5}(4:end))];
    end
    unique_datasets = unique(dataset_nums);
    num_params = max(param_indices);
    
    fprintf('Found data for %d unique datasets\n', length(unique_datasets));
    fprintf('Number of parameter combinations: %d\n\n', num_params);
    
    % Initialize arrays for metrics
    num_datasets = length(unique_datasets);
    metrics.kernel_quality = nan(num_datasets, num_params);
    metrics.activation_accuracy = nan(num_datasets, num_params);
    metrics.runtime = nan(num_datasets, num_params);
    metrics.residuals = cell(num_datasets, num_params);
    metrics.relative_changes = cell(num_datasets, num_params);
    metrics.noise_levels = nan(num_datasets, 1);
    metrics.param_combinations = [];
    
    % Initialize progress bar
    fprintf('Loading and processing files:\n');
    progress_str = '';
    
    % Load each file and extract metrics
    for i = 1:length(files)
        % Update progress bar
        fprintf(repmat('\b', 1, length(progress_str)));
        progress_str = sprintf('Processing file %d/%d...', i, length(files));
        fprintf(progress_str);
        
        % Load result file
        data = load(fullfile(folder_path, files(i).name));
        
        % Extract dataset number and parameter index
        parts = split(files(i).name, {'_', '.'});
        dataset_idx = find(unique_datasets == str2double(parts{3}(8:end)));
        param_idx = str2double(parts{5}(4:end));
        
        % Store parameter combinations if not already stored
        if isempty(metrics.param_combinations)
            metrics.param_combinations = data.param_combinations;
        end
        
        % Extract metrics from extras
        extras = data.extras{1};  % Unwrap from cell
        
        % Store kernel quality (final value)
        if isfield(extras.phase1, 'kernel_quality_factors')
            kq = extras.phase1.kernel_quality_factors;
            metrics.kernel_quality(dataset_idx, param_idx) = kq(end);
        end
        
        % Store activation accuracy (final value)
        if isfield(extras.phase1, 'activation_metrics')
            aa = extras.phase1.activation_metrics;
            metrics.activation_accuracy(dataset_idx, param_idx) = aa(end);
        end
        
        % Store runtime
        if isfield(extras, 'runtime')
            metrics.runtime(dataset_idx, param_idx) = extras.runtime;
        end
        
        % Calculate and store residuals
        if isfield(extras, 'residuals')
            metrics.residuals{dataset_idx, param_idx} = extras.residuals;
        end
        
        % Calculate relative changes
        if isfield(extras, 'residuals')
            res = extras.residuals;
            rel_changes = abs(diff(res)) ./ abs(res(1:end-1));
            metrics.relative_changes{dataset_idx, param_idx} = rel_changes;
        end
        
        % Store noise level (assuming it's in extras)
        if isfield(extras, 'noise_level') && isnan(metrics.noise_levels(dataset_idx))
            metrics.noise_levels(dataset_idx) = extras.noise_level;
        end
    end
    
    % Clear progress bar and print completion
    fprintf('\nData loading complete!\n\n');
    
    % Add metadata
    metrics.num_datasets = num_datasets;
    metrics.num_params = num_params;
    metrics.lambda1_values = unique(metrics.param_combinations(:,1));
    metrics.mini_loop_values = unique(metrics.param_combinations(:,2));
    
    % Print summary
    fprintf('Summary:\n');
    fprintf('- Number of datasets: %d\n', metrics.num_datasets);
    fprintf('- Lambda1 values: [%s]\n', join(string(metrics.lambda1_values), ', '));
    fprintf('- Mini-loop values: [%s]\n', join(string(metrics.mini_loop_values), ', '));
end