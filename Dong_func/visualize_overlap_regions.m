function visualize_overlap_regions(Xiter, k)
    % Get dimensions
    num_kernels = size(Xiter, 3);
    
    % Calculate overlap regions
    IDX = calculate_overlap_regions(Xiter, k);
    
    % Calculate subplot layout
    num_cols = ceil(sqrt(num_kernels));
    num_rows = ceil(num_kernels/num_cols);
    
    % Visualize
    figure('Position', [100, 100, 200*num_cols, 200*num_rows*2]);
    
    % Plot activation regions
    for n = 1:num_kernels
        subplot(2*num_rows, num_cols, n);
        
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
        subplot(2*num_rows, num_cols, n + num_cols*num_rows);
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