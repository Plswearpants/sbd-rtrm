function [Xsol, info] = Xsolve_FISTA(Y, A, lambda, mu, varargin)
%XSOLVE_FISTA   Solve for X using FISTA method
%   - Core usage:
%       [Xsol, info] = Xsolve_FISTA(Y, A, lambda, mu)
%
%   - Optional variables:
% %       [...] = Xsolve_FISTA(..., Xinit, Xpos, getbias)
%       Xinit:      initial value for X (cell array for multiple kernels)
%       Xpos:       constrain X to be a positive solution
%       getbias:    extract constant bias as well as X

    % Initialize variables and function handles:
    fpath = fileparts(mfilename('fullpath'));
    addpath([fpath '/helpers']);
    load([fpath '/../config/Xsolve_config.mat']); %#ok<*LOAD>
    g = huber(mu);

    m = size(Y);
    if (numel(m) > 2)
        n = m(3); m = m(1:2);
    else
        n = 1;
    end

    % Ensure A is always a cell array
    if ~iscell(A)
        A = {A};
    end
    num_kernels = length(A);
    %% Checking arguments:
    nvararg = numel(varargin);
    if nvararg > 3
        error('Too many input arguments.');
    end

    idx = 1; X = cell(1, num_kernels); b = zeros(n, num_kernels);
    if nvararg >= idx && ~isempty(varargin{idx})
        X = varargin{idx}.X;
        b = varargin{idx}.b;
    else
        for i = 1:num_kernels
            X{i} = zeros(m);
        end
    end

    idx = 2; xpos = false;
    if nvararg >= idx && ~isempty(varargin{idx})
        xpos = varargin{idx};
    end
    
    idx = 3; getbias = false;
    if nvararg >= idx && ~isempty(varargin{idx})
        getbias = varargin{idx};
    end

    %% Iterate:    
    t = 1; W = X; u = b;
    costs = NaN(MAXIT, 2);
    doagain = true; it = 0; count = 0;
    while doagain
        it = it + 1;
        % Gradients and Hessians:
        grad_fW = cell(1, num_kernels);
        grad_fu = zeros(n, num_kernels);
        R_A = cell(1, num_kernels);
        for k = 1:num_kernels
            grad_fW{k} = zeros(m);
            R_A{k} = zeros(m);
            for i = 1:n     % sum up
                % Use the entire A{k} matrix instead of indexing it
                Ri = convfft2(A{k}, W{k}(:,:,i)) + u(i,k) - Y(:,:,i);
                grad_fW{k} = grad_fW{k} + convfft2(A{k}, Ri, 1);
                grad_fu(i,k) = sum(Ri(:));
                R_A{k} = R_A{k} + abs(fft2(A{k}, m(1), m(2))).^2;
            end
        end

        % FISTA update
        X_ = cell(1, num_kernels);
        for k = 1:num_kernels
            L = max(R_A{k}(:));
            
            % Debug information
            fprintf('Kernel %d:\n', k);
            fprintf('Size of W{k}: %s\n', mat2str(size(W{k})));
            fprintf('Size of grad_fW{k}: %s\n', mat2str(size(grad_fW{k})));
            
            % Ensure W{k} and grad_fW{k} have the same dimensions
            if ~isequal(size(W{k}), size(grad_fW{k}))
                error('Dimension mismatch between W{k} and grad_fW{k} for kernel %d', k);
            end
            
            % Remove squeeze if W{k} is already 2D
            if ndims(W{k}) == 2
                X_{k} = g.prox(W{k} - 1/L*grad_fW{k}, lambda/L, xpos);
            else
                X_{k} = g.prox(squeeze(W{k} - 1/L*grad_fW{k}), lambda/L, xpos);
            end
            
            % Debug information
            fprintf('Size of X_{k}: %s\n', mat2str(size(X_{k})));
        end
        t_ = (1+sqrt(1+4*t^2))/2;
        for k = 1:num_kernels
            W{k} = X_{k} + (t-1)/t_*(X_{k}-X{k});
        end
        if getbias
            b_ = u - grad_fu/(2*prod(m)*sqrt(n));  
            u = b_ + (t-1)/t_*(b_-b);
            b = b_;
        end
        X = X_; t = t_;

        % Check conditions to repeat iteration:
        f = 0;
        for i = 1:n
            sum_conv = zeros(m);
            for k = 1:num_kernels
                sum_conv = sum_conv + convfft2(A{k}(:,:,i), reshape(X{k}, m)) + b(i,k);
            end
            f = f + norm(sum_conv - Y(:,:,i), 'fro')^2/2;
        end
        costs(it,1) = f;
        costs(it,2) = sum(cellfun(@(x) g.cost(x, lambda), X));

        delta = 0;
        for k = 1:num_kernels
            tmp = grad_fW{k};
            for i = 1:n
                tmp = tmp + grad_fu(i,k);
            end
            delta = delta + norm(g.diffsubg(X{k}, -tmp, lambda, xpos), 'fro')^2;
        end
        delta = sqrt(delta)/sqrt(prod(m)*num_kernels);
        
        if delta < EPSILON
            count = count + 1;
        else
            count = 0;
        end
        doagain = count < 10 && (it < MAXIT);
    end

    % Ensure output is consistent with input
    if num_kernels == 1
        Xsol.X = Xsol.X{1};
        Xsol.W = Xsol.W{1};
        Xsol.b = Xsol.b(:,1);
    end

    % Return solution:
    Xsol.X = X;
    Xsol.b = b;
    Xsol.W = W;
    Xsol.f = sum(costs(it,:));
    info.numit = it;
    info.delta = delta;
    info.costs = costs(1:it,:);
end
