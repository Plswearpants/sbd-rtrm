function [activation_similarity, kernel_similarity] = computeQualityMetrics(X0, Xout, A0, Aout, kernel_size)
    % Bundles both activation and kernel quality metrics into one function
    % without visualization
    %
    % Inputs:
    %   X0: ground truth activation maps [height x width x num_kernels]
    %   Xout: reconstructed activation maps [height x width x num_kernels]
    %   A0: ground truth kernels (cell array)
    %   Aout: reconstructed kernels (cell array)
    %   kernel_size: size of each kernel [num_kernels x 2]
    %
    % Outputs:
    %   activation_similarity: array of similarity scores for activations
    %   kernel_similarity: array of similarity scores for kernels

    % Compute activation metrics (with visualization off)
    activation_metric = evaluateActivationReconstruction(X0, Xout, kernel_size, false);
    activation_similarity = activation_metric.similarity;
    
    % Compute kernel metrics (with visualization off)
    kernel_metric = evaluateKernelQuality(Aout, A0, false);
    kernel_similarity = kernel_metric;
end 