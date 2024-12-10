function [demixing_score, correlation_matrix] = computeDemixingMetric(X)
    % Compute demixing metric for activation maps
    %
    % Input:
    %   X: activation maps [height x width x num_kernels]
    %
    % Output:
    %   demixing_score: scalar value representing overall separation quality
    %                   (higher is better - means better separation)
    %   correlation_matrix: matrix of cross-correlations between activations
    %                      (excluding self-correlations)
    
    num_kernels = size(X, 3);
    correlation_matrix = zeros(num_kernels, num_kernels);
    
    % Compute cross-correlations between all pairs
    for i = 1:num_kernels
        X1 = X(:,:,i);
        for j = i+1:num_kernels
            X2 = X(:,:,j);
            
            % Compute normalized cross-correlation
            xcorr = normxcorr2(X1, X2);
            
            % Take maximum absolute correlation value
            max_corr = max(abs(xcorr(:)));
            
            % Store in symmetric matrix
            correlation_matrix(i,j) = max_corr;
            correlation_matrix(j,i) = max_corr;
        end
    end
    
    % Calculate overall demixing score
    % First get average correlation
    avg_correlation = sum(correlation_matrix(:)) / (num_kernels * (num_kernels - 1));
    % Convert to separation score (1 = perfect separation, 0 = complete correlation)
    demixing_score = 1 - avg_correlation;
end 