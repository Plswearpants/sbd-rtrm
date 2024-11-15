function [ Aout, Xout, bout, extras ] = SBD_test_multi_demixing( Y, k, params, dispfun, kernel_initialguess, AX_iteration, maxIT)
    %SBD Summary of this function goes here
    %
    %   PARAMS STRUCT:
    %   ===============
    %   The options struct should include the fields:
    %       lambda1,  float > 0  : regularization parameter for Phase I
    %       phase2,   bool       : whether to do Phase II (refinement) or not
    %
    %   IF phase2 == true, then the following fields should also be included:
    %       kplus,    int > 0    : border padding (pixels) for sphere lifting
    %       lambda2,  float > 0  : FINAL reg. param. value for Phase II
    %
    %       nrefine,  int >= 1   : number of refinements for Phase II.
    %           Refinement 1 lifts the sphere and uses lambda1, successive
    %           refinements decrease lambda down to lambda2;
    %           i.e. if nrefine == 1, then no decrease in lambda is made.
    %
    %   Finally, two optional fields for the struct. These features are
    %   automatically disabled if the fields are not included or are empty:
    %
    %       Xsolve,   string     : Pick which Xsolve to use--'FISTA' or 'pdNCG.'
    %       xpos,     bool       : Constrain X to have nonnegative entries
    %       getbias,  bool       : Extract constant bias from observation.
    %                                                                                                                      
    %   k: n*2 matrix containing a list of n different kernel sizes [x,y]. 
    %   kernel_initialguess: cell array containing different intial guess of kernel 
    
    %% Process input arguments
    if nargin < 4 || isempty(dispfun)
        dispfun = @(Y,A,X,k,kplus,idx) 0;
    end
    
    lambda1 = params.lambda1;
    if params.phase2
        kplus = params.kplus;
        lambda2 = params.lambda2;
        nrefine = params.nrefine;
    end
    
    if ~isfield(params, 'xpos') || isempty(params.xpos)
        xpos = false;
    else
        xpos = params.xpos;
    end
    
    if ~isfield(params, 'getbias') || isempty(params.getbias)
        getbias = false;
    else
        getbias = params.getbias;
    end
    
    if ~isfield(params, 'Xsolve') || isempty(params.Xsolve)
        Xsolve = 'FISTA';
    else
        Xsolve = params.Xsolve;
    end
    
    kernel_num = size(k,1);
    mu = 10^-6;
    
    % Update the configuration file with the new max_iteration
    update_config('Xsolve_config.mat', 'MAXIT', AX_iteration, 'Xsolve_config_tunable.mat');
    update_config('Asolve_config.mat','options.maxiter', AX_iteration, 'Asolve_config_tunable.mat');
    
    % Extract X0 and A0 from params
    if ~isfield(params, 'X0') || ~isfield(params, 'A0')
        error('params must contain X0 and A0 for quality metrics');
    end
    X0 = params.X0;
    A0 = params.A0;
    
    %% Phase I: Initialization and First Iteration
    fprintf('PHASE I: Initialization and First Iteration\n');
    A = kernel_initialguess;
    X_struct = struct();
    Xiter = zeros([size(Y),kernel_num]);
    biter = zeros(kernel_num,1);
    
    % Initialize quality metrics arrays in extras
    extras.phase1.activation_metrics = cell(1, maxIT);
    extras.phase1.kernel_quality_factors = zeros(maxIT, kernel_num);
    
    % Add demixing factor (only for Phase I)
    faint_factor = 1;
    
    for iter = 1:maxIT
        starttime = tic;
    
        % Compute Y_background using demixing approach
        Y_sum = zeros(size(Y));
        for m = 1:kernel_num
            if iter > 1
                Y_sum = Y_sum + convfft2(A{m}, Xiter(:,:,m));
            end
        end
        Y_residual = Y - Y_sum;
        
        % Update each kernel with demixing
        for n = 1:kernel_num
            % Calculate Yiter for this kernel using demixing approach
            Yiter = Y_residual + (iter > 1) * (1-1/(faint_factor*iter+1))*convfft2(A{n}, Xiter(:,:,n)) + (1/(faint_factor*iter+1))*Y_sum;
            
            dispfun1 = @(A, X) dispfun{n}(Y, A, X, k(n,:), []);
            
            if iter == 1
                % Initial X computation
                X_struct.(['x',num2str(n)]) = Xsolve_FISTA_tunable(Y, A{n}, lambda1(n), mu, [], xpos);
            end
            
            [A{n}, X_struct.(['x',num2str(n)]), info] = Asolve_Manopt_tunable(Yiter, A{n}, lambda1(n), Xsolve, X_struct.(['x',num2str(n)]), xpos, getbias, dispfun1);
            
            Xiter(:,:,n) = X_struct.(['x',num2str(n)]).X;
            biter(n) = X_struct.(['x',num2str(n)]).b;
        end
    
        % Keep the same quality metrics computation
        [activation_metrics, kernel_metrics] = computeQualityMetrics(X0, Xiter, A0, A, k);
        extras.phase1.activation_metrics{iter} = activation_metrics;
        extras.phase1.kernel_quality_factors(iter,:) = kernel_metrics;
        
        runtime = toc(starttime);
        fprintf('Iteration %d: Runtime = %.2fs\n', iter, runtime);
    end
    
    % Store results
    extras.phase1.Aout = A;
    extras.phase1.Xout = Xiter;
    extras.phase1.biter = biter;
    
    %% PHASE II: Refinement (if requested) - Using original approach without demixing
    if params.phase2
        fprintf('\nPHASE II: Refinement\n');
        k3 = k + 2*kplus;
    
        A2 = cell(1, kernel_num);
        X2_struct = struct();
        for n = 1:kernel_num
            A2{n} = zeros(k3(n,:));
            A2{n}(kplus(n,1)+(1:k(n,1)), kplus(n,2)+(1:k(n,2))) = A{n};
            X2_struct.(['x',num2str(n)]) = X_struct.(['x',num2str(n)]);
        end
    
        lambda = lambda1;
        lam2fac = (lambda2./lambda1).^(1/nrefine);
        
        % Initialize Phase II metrics
        extras.phase2.activation_metrics = cell(1, nrefine + 1);
        extras.phase2.kernel_quality_factors = zeros(nrefine + 1, kernel_num);
        
        for i = 1:nrefine + 1
            fprintf('lambda iteration %d/%d: \n', i, nrefine + 1);
            
            % Standard Y_residual calculation (no demixing)
            for m = 1:kernel_num
                Y_sum = Y_sum + convfft2(A2{m}, X2_struct.(['x',num2str(m)]).X);
            end
            Y_residual = Y - Y_sum;
    
            for n = 1:kernel_num
                fprintf('Processing kernel %d, lambda = %.1e: \n', n, lambda(n));
                % Calculate Yiter without demixing
                Yiter = Y_residual + convfft2(A2{n}, X2_struct.(['x',num2str(n)]).X);
                
                dispfun2 = @(A, X) dispfun{n}(Y, A, X, k3(n,:), kplus(n,:));
                [A2{n}, X2_struct.(['x',num2str(n)]), info] = Asolve_Manopt_tunable(Yiter, A2{n}, lambda(n), Xsolve, X2_struct.(['x',num2str(n)]), xpos, getbias, dispfun2);
                
                % Attempt to 'unshift" the a and x
                score = zeros(2*kplus(n,1)+1, 2*kplus(n,2)+1);
                for tau1 = -kplus(n,1):kplus(n,1)
                    ind1 = tau1+kplus(n,1)+1;
                    for tau2 = -kplus(n,2):kplus(n,2)
                        ind2 = tau2+kplus(n,2)+1;
                        temp = A2{n}(ind1:(ind1+k(n,1)-1), ind2:(ind2+k(n,2)-1));
                        score(ind1,ind2) = norm(temp(:), 1);
                    end
                end
                [temp,ind1] = max(score); [~,ind2] = max(temp);
                tau = [ind1(ind2) ind2]-kplus(n,:)-1;
                A2{n} = circshift(A2{n},-tau);
                X2_struct.(['x',num2str(n)]).X = circshift(X2_struct.(['x',num2str(n)]).X,tau);
                X2_struct.(['x',num2str(n)]).W = circshift(X2_struct.(['x',num2str(n)]).W,tau);
    
                % Save phase 2 extras
                extras.phase2.A{n} = A2{n};
                extras.phase2.X{n} = X2_struct.(['x',num2str(n)]).X;
                extras.phase2.b(n) = X2_struct.(['x',num2str(n)]).b;
                extras.phase2.info{n} = info;
            end
            
            % Evaluate metrics for this refinement
            X2_combined = zeros(size(Y,1), size(Y,2), kernel_num);
            A2_central = cell(1, kernel_num);
            for n = 1:kernel_num
                X2_combined(:,:,n) = X2_struct.(['x',num2str(n)]).X;
                A2_central{n} = A2{n}(kplus(n,1)+(1:k(n,1)), kplus(n,2)+(1:k(n,2)));
            end
    
            [activation_metrics, kernel_metrics] = computeQualityMetrics(X0, X2_combined, A0, A2_central, k3);
            extras.phase2.activation_metrics{i} = activation_metrics;
            extras.phase2.kernel_quality_factors(i,:) = kernel_metrics;
            
            % Print results
            fprintf('Refinement %d Metrics:\n', i);
            fprintf('Activation Quality:\n');
            for n = 1:kernel_num
                fprintf('Kernel %d - Similarity: %.3f\n', n, activation_metrics(n));
            end
            fprintf('Kernel Quality Factors:\n');
            for n = 1:kernel_num
                fprintf('Kernel %d: %.3f\n', n, kernel_metrics(n));
            end
            
            lambda = lambda .* lam2fac;
        end
    end
    
    %% Finished: get the final A, X
    if params.phase2
        Aout = cell(1, kernel_num);
        Xout = zeros(size(Y,1), size(Y,2), kernel_num);
        bout = zeros(kernel_num, 1);
        extras.normA = zeros(kernel_num, 1);
        
        for n = 1:kernel_num
            Aout{n} = A2{n}(kplus(n,1)+(1:k(n,1)), kplus(n,2)+(1:k(n,2)));
            extras.normA(n) = norm(Aout{n}(:));
            Xout(:,:,n) = circshift(X2_struct.(['x',num2str(n)]).X, kplus(n,:)) * extras.normA(n);
            Aout{n} = Aout{n} / extras.normA(n);
            bout(n) = X2_struct.(['x',num2str(n)]).b;
        end
    else
        Aout = A;
        extras.normA = zeros(kernel_num, 1);
        for n = 1:kernel_num
            extras.normA(n) = norm(Aout{n}(:));
        end
        Xout = Xiter;
        bout = biter;
    end
    
    runtime = toc(starttime);
    fprintf('\nDone! Runtime = %.2fs. \n\n', runtime);
    
    % Call the visualization function
    visualize_quality_factors(extras, params.phase2, maxIT, nrefine, kernel_num);
    end
    
    function visualize_quality_factors(extras, phase2_performed, maxIT, nrefine, kernel_num)
        figure;
    
        % Plot Activation Similarity
        subplot(2,1,1);
        similarities_phase1 = zeros(maxIT, kernel_num);
        for i = 1:maxIT
            similarities_phase1(i,:) = extras.phase1.activation_metrics{i};
        end
        
        colors = lines(kernel_num);
        for n = 1:kernel_num
            plot(1:maxIT, similarities_phase1(:,n), 'Color', colors(n,:), 'LineStyle', '-', ...
                'DisplayName', sprintf('Phase I Kernel %d', n));
            hold on;
        end
        
        if phase2_performed
            similarities_phase2 = zeros(nrefine+1, kernel_num);
            for i = 1:nrefine+1
                similarities_phase2(i,:) = extras.phase2.activation_metrics{i};
            end
            
            for n = 1:kernel_num
                plot(maxIT:maxIT+nrefine, similarities_phase2(:,n), 'Color', colors(n,:), ...
                    'LineStyle', ':', 'DisplayName', sprintf('Phase II Kernel %d', n));
            end
        end
        xlabel('Iteration');
        ylabel('Activation Similarity');
        title('Activation Reconstruction Quality vs. Iteration');
        legend('show');
        grid on;
        hold off;
    
        % Plot Kernel Quality Factors
        subplot(2,1,2);
        for n = 1:kernel_num
            plot(1:maxIT, extras.phase1.kernel_quality_factors(:,n), 'Color', colors(n,:), ...
                'LineStyle', '-', 'DisplayName', sprintf('Phase I Kernel %d', n));
            hold on;
        end
        
        if phase2_performed
            for n = 1:kernel_num
                plot(maxIT:maxIT+nrefine, extras.phase2.kernel_quality_factors(:,n), ...
                    'Color', colors(n,:), 'LineStyle', ':', ...
                    'DisplayName', sprintf('Phase II Kernel %d', n));
            end
        end
        xlabel('Iteration');
        ylabel('Kernel Quality Factor');
        title('Kernel Quality Factors vs. Iteration');
        legend('show');
        grid on;
        hold off;
    end