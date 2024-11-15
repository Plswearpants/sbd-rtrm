function [metrics, aligned_maps] = evaluateActivationReconstruction(X0, Xout, kernel_size, visualize)
    % Evaluate reconstruction quality by first aligning activation maps and then
    % computing similarity with density-adaptive Gaussian filtering
    %
    % Inputs:
    %   X0: original activation maps [height x width x num_kernels]
    %   Xout: reconstructed activation maps [height x width x num_kernels]
    %   kernel_size: size of each kernel [num_kernels x 2]
    %   visualize: (optional) boolean for visualization (must be 1 or 0)
    %
    % Outputs:
    %   metrics: similarity scores for each kernel (if visualize=false)
    %           or full metrics struct (if visualize=true)
    %   aligned_maps: aligned activation maps and filtered results
    
    if nargin < 4
        visualize = false;
    end
    
    % Step 1: Align activation maps
    [Xout_aligned, offsets, align_quality] = alignActivationMaps(X0, Xout, kernel_size);
    
    % Step 2: Compute similarity with density-adaptive Gaussian
    [similarities, filtered_maps] = computeActivationSimilarity(X0, Xout_aligned, kernel_size, visualize);
    
    % Store aligned maps
    aligned_maps.Xout_aligned = Xout_aligned;
    aligned_maps.filtered = filtered_maps;
    
    if ~visualize
        metrics.similarity = similarities;
    else
        % Full metrics output and visualization when requested
        metrics.alignment = align_quality;
        metrics.offset = offsets;
        metrics.similarity = similarities;
        metrics.density = [filtered_maps.density];
        metrics.L = [filtered_maps.L];
        metrics.sigma = [filtered_maps.sigma];
        
        % Print results only if visualization is true
        fprintf('\nFinal Evaluation Results:\n');
        fprintf('========================\n');
        for k = 1:num_kernels
            fprintf('\nKernel %d:\n', k);
            fprintf('  Pattern Match: %.3f\n', metrics(k).alignment.primary_peak);
            fprintf('  Similarity: %.3f\n', metrics(k).similarity);
        end
    end
end 