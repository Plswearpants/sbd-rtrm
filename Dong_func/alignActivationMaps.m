function [Xout_aligned, offset, quality] = alignActivationMaps(X0, Xout, kernel_size)
    % Input:
    %   X0: original activation maps [height x width x num_kernels]
    %   Xout: reconstructed activation maps [height x width x num_kernels]
    %   kernel_size: size of each kernel [num_kernels x 2]
    
    num_kernels = size(X0, 3);
    Xout_aligned = zeros(size(Xout));
    offset = zeros(num_kernels, 2);
    
    % Create figure for visualizations
    fig = figure('Name', 'Activation Maps Alignment');
    
    for k = 1:num_kernels
        fprintf('\nProcessing Kernel %d:\n', k);
        
        % Broaden activation maps with Gaussian window
        sigma = min(kernel_size(k,:)) / 10;
        window_size = ceil(3 * sigma);
        [x, y] = meshgrid(-window_size:window_size);
        gaussian_kernel = exp(-(x.^2 + y.^2)/(2*sigma^2));
        gaussian_kernel = gaussian_kernel / sum(gaussian_kernel(:));
        
        % Apply Gaussian broadening
        X0_broad = conv2(X0(:,:,k), gaussian_kernel, 'same');
        Xout_broad = conv2(Xout(:,:,k), gaussian_kernel, 'same');
        
        % Normalize broadened maps
        X0_broad = X0_broad / max(X0_broad(:));
        Xout_broad = Xout_broad / max(Xout_broad(:));
        
        % Apply normalized cross-correlation
        c = normxcorr2(X0_broad, Xout_broad);
        
        % Find peak correlation and secondary peaks
        [max_corr, max_idx] = max(c(:));
        [y_peak, x_peak] = ind2sub(size(c), max_idx);
        
        % Calculate offset
        y_offset = y_peak - size(X0, 1);
        x_offset = x_peak - size(X0, 2);
        
        % Store offset
        offset(k,:) = [y_offset, x_offset];
        
        % Apply shift to original Xout
        Xout_aligned(:,:,k) = circshift(Xout(:,:,k), [-y_offset, -x_offset]);
        
        % Compute quality metrics for this kernel
        quality(k).primary_peak = max_corr;
        quality(k).offset = [y_offset, x_offset];
        
        % Find secondary peaks
        c_copy = c;
        peak_region = ceil(min(kernel_size(k,:))/4);
        y_range = max(1, y_peak-peak_region):min(size(c,1), y_peak+peak_region);
        x_range = max(1, x_peak-peak_region):min(size(c,2), x_peak+peak_region);
        c_copy(y_range, x_range) = -inf;
        [secondary_peak, ~] = max(c_copy(:));
        
        quality(k).peak_ratio = secondary_peak / max_corr;
        quality(k).mean_correlation = mean(c(:));
        quality(k).peak_sharpness = max_corr / quality(k).mean_correlation;
        
        % Visualize initial, aligned, and original activation maps
        subplot(num_kernels, 3, 3*(k-1) + 1);
        imagesc(Xout(:,:,k));
        title(sprintf('K%d: Initial Xout', k));
        colorbar;
        axis square;
        
        subplot(num_kernels, 3, 3*(k-1) + 2);
        imagesc(Xout_aligned(:,:,k));
        title(sprintf('K%d: Aligned Xout\nMatch: %.2f | Unique: %.2f | Precise: %.2f', ...
            k, quality(k).primary_peak, ...    % How well patterns match
            1 - quality(k).peak_ratio, ...     % Inverted for intuition (1 = unique)
            quality(k).peak_sharpness));       % How precise the alignment is
        colorbar;
        axis square;
        
        subplot(num_kernels, 3, 3*(k-1) + 3);
        imagesc(X0(:,:,k));
        title(sprintf('K%d: Original X0', k));
        colorbar;
        axis square;
        
        % Print quality metrics with interpretations
        fprintf('Kernel %d Quality Metrics:\n', k);
        fprintf('  Pattern Match: %.3f  (1 = perfect match)\n', quality(k).primary_peak);
        fprintf('  Uniqueness: %.3f     (1 = single solution)\n', 1 - quality(k).peak_ratio);
        fprintf('  Precision: %.3f      (>1 = precise alignment)\n', quality(k).peak_sharpness);
        fprintf('  Offset: [%d, %d]\n\n', offset(k,1), offset(k,2));
    end
    
    % Adjust figure layout
    set(fig, 'Position', [100, 100, 1200, 300*num_kernels]);
    sgtitle('Activation Maps Alignment (Initial → Aligned → Original)');
end
