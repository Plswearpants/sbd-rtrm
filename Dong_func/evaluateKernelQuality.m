function quality_factors = evaluateKernelQuality(output_kernels, gt_kernels, do_visualization)
    % Evaluates kernel quality by comparing QPI patterns
    %
    % Inputs:
    %   output_kernels: Single kernel or cell array of kernels from algorithm
    %   gt_kernels: Single kernel or cell array of ground truth kernels
    %   do_visualization: Boolean flag for visualization (default: false)
    %
    % Outputs:
    %   quality_factors: Array of similarity scores for each kernel pair
    
    % Handle default visualization parameter
    if nargin < 3
        do_visualization = false;
    end
    
    % Convert to cell if single kernel
    if ~iscell(output_kernels)
        output_kernels = {output_kernels};
    end
    if ~iscell(gt_kernels)
        gt_kernels = {gt_kernels};
    end
    
    num_kernels = length(output_kernels);
    assert(length(gt_kernels) == num_kernels, 'Number of output and GT kernels must match');
    
    % Initialize quality factors array
    quality_factors = zeros(1, num_kernels);
    
    if do_visualization
        % Create figure with appropriate subplots
        figure('Position', [100, 100, 1200, 300*num_kernels]);
    end
    
    for k = 1:num_kernels
        [score, qpi_out, qpi_gt] = kernel_QPI_metric(output_kernels{k}, gt_kernels{k});
        quality_factors(k) = score;
        
        if do_visualization
            % Plot original kernels
            subplot(num_kernels, 4, (k-1)*4 + 1);
            imagesc(output_kernels{k});
            title(sprintf('Output Kernel %d', k));
            axis square; colorbar;
            
            subplot(num_kernels, 4, (k-1)*4 + 2);
            imagesc(gt_kernels{k});
            title(sprintf('GT Kernel %d', k));
            axis square; colorbar;
            
            % Plot QPI patterns
            subplot(num_kernels, 4, (k-1)*4 + 3);
            imagesc(log(qpi_out + 1));  % log scale for better visualization
            title(sprintf('Output QPI %d (Score: %.3f)', k, score));
            axis square; colorbar;
            
            subplot(num_kernels, 4, (k-1)*4 + 4);
            imagesc(log(qpi_gt + 1));   % log scale for better visualization
            title(sprintf('GT QPI %d', k));
            axis square; colorbar;
        end
    end
    
    if do_visualization && num_kernels > 1
        mean_score = mean(quality_factors);
        sgtitle(sprintf('Kernel QPI Comparison (Mean Score: %.3f)', mean_score));
    elseif do_visualization
        sgtitle(sprintf('Kernel QPI Comparison (Score: %.3f)', quality_factors(1)));
    end
end 