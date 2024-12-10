% Load the metrics
metrics = load_parallel_results();

%% Create the heatmaps
plot_performance_heatmaps(metrics);

%% Create parameter slices  
plot_parameter_slices(metrics);

%% Single dataset visualization
plot_single_dataset(metrics);

%% Optional: Save the figures
%saveas(gcf, 'parameter_slices.fig');
%saveas(gcf, 'parameter_slices.png'); 