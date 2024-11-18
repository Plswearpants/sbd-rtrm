function visualizeResults(Y, A0, Aout, X0, Xout, bout, extras)
    % Visualize results from SBD-STM reconstruction
    %
    % Inputs:
    %   Y: Original observation
    %   A0: Ground truth kernels (cell array)
    %   Aout: Reconstructed kernels (cell array)
    %   X0: Ground truth activation maps
    %   Xout: Reconstructed activation maps
    %   bout: Bias terms [num_kernels x 1]
    %   extras: Struct containing convergence metrics
    
    num_kernels = length(A0);
    
    % Get kernel sizes from A0
    kernel_size = zeros(num_kernels, 2);
    for n = 1:num_kernels
        kernel_size(n,:) = size(A0{n});
    end
    
    % Compute reconstructed image
    Y_reconstructed = zeros(size(Y));
    for i = 1:num_kernels
        Y_reconstructed = Y_reconstructed + convfft2(Aout{i}, Xout(:,:,i));
        Y_reconstructed = Y_reconstructed + bout(i);  % Add bias term
    end

    % 1. Original vs Reconstructed Full Image
    figure('Name', 'Full Image Comparison');
    subplot(121);
    imagesc(Y);
    title('Original Image');
    colorbar;
    axis image;

    subplot(122);
    imagesc(Y_reconstructed);
    title('Reconstructed Image');
    colorbar;
    axis image;
    sgtitle('Original vs Reconstructed Image');

    % 2. Kernel Quality Analysis (using QPI patterns)
    fprintf('\nAnalyzing Kernel Quality with QPI Patterns:\n');
    quality_factors = evaluateKernelQuality(Aout, A0, true);
    fprintf('QPI-based Quality Factors:\n');
    for n = 1:num_kernels
        fprintf('Kernel %d: %.4f\n', n, quality_factors(n));
    end

    % 3. Activation Map Analysis
    fprintf('\nAnalyzing Activation Map Quality:\n');
    [activation_metrics, aligned_maps] = evaluateActivationReconstruction(X0, Xout, kernel_size, true);

    % Print detailed activation metrics
    fprintf('\nActivation Reconstruction Metrics:\n');
    for k = 1:num_kernels
        fprintf('Kernel %d:\n', k);
        fprintf('  Alignment Quality: %.4f\n', activation_metrics.alignment(k).primary_peak);
        fprintf('  Similarity Score: %.4f\n', activation_metrics.similarity(k));
        fprintf('  Offset: [%d, %d]\n', activation_metrics.offset(k,1), activation_metrics.offset(k,2));
    end

    % 4. Convergence History
    figure('Name', 'Convergence History');
    
    % Check if Phase II was performed
    phase2_performed = isfield(extras, 'phase2');
    maxIT = size(extras.phase1.kernel_quality_factors, 1);
    
    % Plot Activation Metrics
    subplot(2,1,1);
    % Phase I
    plot(1:maxIT, extras.phase1.activation_metrics, 'LineStyle', '-');
    hold on;
    % Phase II (if performed)
    if phase2_performed
        nrefine = size(extras.phase2.activation_metrics, 1) - 1;  % Subtract 1 because nrefine+1 is total iterations
        plot(maxIT:(maxIT+nrefine), extras.phase2.activation_metrics, 'LineStyle', ':');
    end
    hold off;
    title('Activation Metrics Convergence');
    xlabel('Iteration');
    ylabel('Activation Similarity');
    legend(arrayfun(@(x) sprintf('Kernel %d', x), 1:num_kernels, 'UniformOutput', false));
    grid on;

    % Plot Kernel Quality Factors
    subplot(2,1,2);
    % Phase I
    plot(1:maxIT, extras.phase1.kernel_quality_factors, 'LineStyle', '-');
    hold on;
    % Phase II (if performed)
    if phase2_performed
        plot(maxIT:(maxIT+nrefine), extras.phase2.kernel_quality_factors, 'LineStyle', ':');
    end
    hold off;
    title('Kernel Quality Factors Convergence');
    xlabel('Iteration');
    ylabel('Quality Factor');
    legend(arrayfun(@(x) sprintf('Kernel %d', x), 1:num_kernels, 'UniformOutput', false));
    grid on;
end 