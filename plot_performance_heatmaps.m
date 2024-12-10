function plot_performance_heatmaps(metrics)
    % Create figure for performance metrics
    fig = figure('Position', [100 100 1400 800], 'Name', 'Performance Metrics');
    
    % Create a panel for controls that will scale with the figure
    control_panel = uipanel('Parent', fig, ...
        'Units', 'normalized', ...
        'Position', [0 0.95 1 0.05], ... % Top 5% of figure
        'BackgroundColor', get(fig, 'Color'));
    
    % Add dataset selector with normalized positions within the panel
    uicontrol('Parent', control_panel, ...
        'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.02 0.1 0.1 0.8], ...
        'String', 'Highlight Dataset:', ...
        'BackgroundColor', get(fig, 'Color'));
    
    dataset_selector = uicontrol('Parent', control_panel, ...
        'Style', 'popup', ...
        'Units', 'normalized', ...
        'Position', [0.13 0.1 0.1 0.8], ...
        'String', arrayfun(@(x) sprintf('Dataset %d', x), ...
        1:metrics.num_datasets, 'UniformOutput', false), ...
        'Callback', @updateHighlight);
    
    % Add transparency slider
    uicontrol('Parent', control_panel, ...
        'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.24 0.1 0.12 0.8], ...
        'String', 'Other Datasets Alpha:', ...
        'BackgroundColor', get(fig, 'Color'));
    
    alpha_slider = uicontrol('Parent', control_panel, ...
        'Style', 'slider', ...
        'Units', 'normalized', ...
        'Position', [0.37 0.1 0.1 0.8], ...
        'Min', 0, 'Max', 1, 'Value', 0.3, ...
        'Callback', @updateHighlight);
    
    % Add description text box after other controls
    description_text = uicontrol('Parent', control_panel, ...
        'Style', 'text', ...
        'Units', 'normalized', ...
        'Position', [0.5 0.1 0.45 0.8], ...  % Right half of control panel
        'String', '', ...
        'BackgroundColor', get(fig, 'Color'), ...
        'HorizontalAlignment', 'left');
    
    % Store description text handle in UserData
    fig.UserData.description_text = description_text;
    
    % Create a panel for plots that will scale with the figure
    plot_panel = uipanel('Parent', fig, ...
        'Units', 'normalized', ...
        'Position', [0 0 1 0.95], ... % Bottom 95% of figure
        'BackgroundColor', get(fig, 'Color'));
    
    % Store plot data in figure
    fig.UserData.metrics = metrics;
    fig.UserData.surfaces = cell(4,1);
    fig.UserData.highlighted_dataset = 1;
    fig.UserData.other_alpha = 0.3;
    
    % Create plots
    titles = {'Kernel Quality', 'Activation Recovery', ...
             'Runtime (s)', 'Demixing Score (higher is better)'};
    fields = {'kernel_quality_final', 'activation_accuracy_final', ...
             'runtime', 'demixing_score'};
    
    % Create subplot layout within the plot panel
    t = tiledlayout(plot_panel, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for i = 1:4
        ax = nexttile(t);  % Explicitly specify the tiledlayout
        fig.UserData.surfaces{i} = plot_metric_surface(ax, metrics.(fields{i}), ...
            metrics.lambda1_values, metrics.mini_loop_values, titles{i}, ...
            fig.UserData.highlighted_dataset, fig.UserData.other_alpha);
        title(ax, titles{i});
    end
    
    sgtitle('Performance Metrics Across Parameter Space', 'FontSize', 14);
    
    % Nested function to update highlighting
    function updateHighlight(~,~)
        highlighted_dataset = dataset_selector.Value;
        other_alpha = alpha_slider.Value;
        
        % Update stored values
        fig.UserData.highlighted_dataset = highlighted_dataset;
        fig.UserData.other_alpha = other_alpha;
        
        % Update description text if available
        if isfield(metrics, 'dataset_descriptions')
            description_text.String = metrics.dataset_descriptions{highlighted_dataset};
        end
        
        % Update each plot
        for plot_idx = 1:4
            data = metrics.(fields{plot_idx});
            surfaces = fig.UserData.surfaces{plot_idx};
            
            % Update surface transparencies
            for surf_idx = 1:length(surfaces)
                if surf_idx == highlighted_dataset
                    surfaces(surf_idx).FaceAlpha = 1;
                else
                    surfaces(surf_idx).FaceAlpha = other_alpha;
                end
            end
        end
    end
end

function surfaces = plot_metric_surface(ax, data, lambda1_values, mini_loop_values, title_str, highlighted_dataset, other_alpha)
    % Create meshgrid for surface
    [X, Y] = meshgrid(lambda1_values, mini_loop_values);
    
    hold on
    
    if isscalar(data)
        % For scalar metrics (like runtime)
        Z = reshape(data, size(X));
        surfaces = surf(X, Y, Z, 'FaceAlpha', 0.7);
    else
        % For array metrics, plot each dataset as a separate surface
        num_datasets = size(data, 1);
        surfaces = gobjects(num_datasets, 1);
        
        % Define distinct colors
        distinct_colors = [
            0.8500    0.3250    0.0980;  % Orange
            0         0.4470    0.7410;  % Blue
            0.9290    0.6940    0.1250;  % Yellow
            0.4940    0.1840    0.5560;  % Purple
            0.4660    0.6740    0.1880;  % Green
            0.6350    0.0780    0.1840;  % Dark Red
            0         0.7500    0.7500;  % Cyan
            0.7500    0         0.7500;  % Magenta
            0.2500    0.2500    0.2500;  % Gray
            0.9500    0.5000    0.2000;  % Light Orange
            0.1000    0.5000    0.5000;  % Teal
            0.5000    0.5000    0;       % Olive
        ];
        
        for i = 1:num_datasets
            dataset_data = data(i,:);
            Z = reshape(dataset_data, size(X));
            
            % Set transparency based on whether this is the highlighted dataset
            if i == highlighted_dataset
                alpha = 1;
            else
                alpha = other_alpha;
            end
            
            surfaces(i) = surf(X, Y, Z, 'FaceAlpha', alpha, ...
                'FaceColor', distinct_colors(i,:), ...
                'EdgeColor', 'k', ...
                'EdgeAlpha', 0.2);
        end
    end
    
    % Format axes
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    set(gca, 'XTick', lambda1_values);
    set(gca, 'YTick', mini_loop_values);
    set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.3f', x), lambda1_values, 'UniformOutput', false));
    set(gca, 'YTickLabel', arrayfun(@(x) sprintf('%d', x), mini_loop_values, 'UniformOutput', false));
    
    xlabel('\lambda_1');
    ylabel('mini\_loop');
    zlabel('Value');
    grid on;
    view(45, 30);
    
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