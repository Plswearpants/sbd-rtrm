function quality_factor = compute_kernel_quality_factors(ground_truth_kernel, input_kernel)
    % Ensure both kernels are column vectors
    ground_truth_kernel = ground_truth_kernel(:);
    input_kernel = input_kernel(:);
    
    % Compute the variance of the difference
    quality_factor = var(ground_truth_kernel - input_kernel);
end

