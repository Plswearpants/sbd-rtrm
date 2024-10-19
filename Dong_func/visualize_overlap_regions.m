function visualize_overlap_regions(Xiter, k)
    % Get dimensions
    num_kernels = size(Xiter, 3);
    
    % Calculate overlap regions
    IDX = calculate_overlap_regions(Xiter, k);
    
    % Visualize
    figure('Position', [100, 100, 1200, 800]);
    
    % Plot activation regions
    for n = 1:num_kernels
        subplot(2, 3, n);
        
        % Create kernel of ones
        A = ones(k(n, 1), k(n, 2));
        
        % Convolve activation map with kernel of ones
        activation_region = conv2(Xiter(:,:,n), A, 'same');
        
        % Normalize for better visualization
        activation_region = activation_region / max(activation_region(:));
        
        imagesc(activation_region);
        title(sprintf('Kernel %d Activation Region', n));
        colorbar;
        axis equal tight;
    end
    
    % Plot overlap regions
    for n = 1:num_kernels
        subplot(2, 3, n+3);
        imagesc(Xiter(:,:,n));
        hold on;
        h = imagesc(IDX(:,:,n));
        set(h, 'AlphaData', 0.3 * IDX(:,:,n));
        colormap(gca, [1 1 1; 1 0 0]);  % White background, red overlay
        title(sprintf('Kernel %d Overlap', n));
        colorbar;
        axis equal tight;
        hold off;
    end
end