function plot_performance_heatmaps(metrics)
    % Create figure for performance metrics
    figure('Position', [100 100 1200 800], 'Name', 'Performance Metrics');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % 1. Kernel Quality Surface
    nexttile
    plot_metric_surface(metrics.kernel_quality_final, metrics.lambda1_values, metrics.mini_loop_values, ...
        'Kernel Quality');
    
    % 2. Activation Recovery Surface
    nexttile
    plot_metric_surface(metrics.activation_accuracy_final, metrics.lambda1_values, metrics.mini_loop_values, ...
        'Activation Recovery');
    
    % 3. Runtime Surface
    nexttile
    plot_metric_surface(metrics.runtime, metrics.lambda1_values, metrics.mini_loop_values, ...
        'Runtime (s)');
    
    % 4. Convergence Success Surface
    nexttile
    convergence_map = calculate_convergence_map(metrics);
    plot_metric_surface(convergence_map, metrics.lambda1_values, metrics.mini_loop_values, ...
        'Convergence Success');
    
    sgtitle('Performance Metrics Across Parameter Space', 'FontSize', 14);
end

function plot_metric_surface(data, lambda1_values, mini_loop_values, title_str)
    % Create meshgrid for surface
    [X, Y] = meshgrid(lambda1_values, mini_loop_values);
    
    % Create surface plots for each dataset
    hold on
    
    if isscalar(data)
        % For scalar metrics (like runtime)
        Z = reshape(data, size(X));
        surf(X, Y, Z, 'FaceAlpha', 0.7);
    else
        % For array metrics, plot each dataset as a separate surface
        num_datasets = size(data, 1);
        colors = winter(num_datasets); % Create colormap for different datasets
        
        for i = 1:num_datasets
            dataset_data = data(i,:);
            Z = reshape(dataset_data, size(X));
            surf(X, Y, Z, 'FaceAlpha', 0.7, 'FaceColor', colors(i,:), 'EdgeColor', 'k');
        end
    end
    
    % Format axes
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    
    % Set ticks
    set(gca, 'XTick', lambda1_values);
    set(gca, 'YTick', mini_loop_values);
    
    % Format tick labels
    set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.3f', x), lambda1_values, 'UniformOutput', false));
    set(gca, 'YTickLabel', arrayfun(@(x) sprintf('%d', x), mini_loop_values, 'UniformOutput', false));
    
    % Labels and title
    xlabel('\lambda_1');
    ylabel('mini\_loop');
    zlabel('Value');
    title(title_str);
    
    % Add grid
    grid on;
    
    % Set view angle
    view(45, 30);
    
    % Add colorbar to show dataset index
    if ~isscalar(data)
        colorbar('Ticks', linspace(0,1,num_datasets), ...
                'TickLabels', arrayfun(@(x) sprintf('Dataset %d', x), 1:num_datasets, 'UniformOutput', false));
    end
    
    hold off
end

function convergence_map = calculate_convergence_map(metrics)
    % Initialize convergence map
    convergence_map = zeros(metrics.num_datasets, metrics.num_params);
    
    % For each dataset and parameter combination
    for i = 1:metrics.num_datasets
        for j = 1:metrics.num_params
            % Get kernel quality trajectory
            kq_traj = metrics.kernel_quality_trajectory{i,j};
            
            if ~isempty(kq_traj)
                % Calculate differences between consecutive values
                diffs = abs(diff(kq_traj));
                
                % Check if differences are monotonically decreasing
                is_monotonic = all(diff(diffs) <= 0);
                
                % TODO: Could be improved by:
                % 1. Considering relative improvement compared to ground truth
                % 2. Adding a threshold for "flatness" at convergence
                % 3. Incorporating activation metrics trajectory
                % 4. Using a more sophisticated convergence criterion
                
                if is_monotonic
                    convergence_map(i,j) = 1;
                else
                    convergence_map(i,j) = 0;
                end
            end
        end
    end
end

function color = get_contrast_color(value, min_val, max_val)
    % Return black or white depending on background intensity
    normalized_value = (value - min_val) / (max_val - min_val);
    if normalized_value > 0.5
        color = 'black';
    else
        color = 'white';
    end
end