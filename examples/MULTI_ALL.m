function [ A_phase_I, Aout, Xout, bout, extras ] = SBD_test( Y, k, params, dispfun, kernel_initialguess )
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
%
%   Finally, two optional fields for the struct. These features are
%   automatically disabled if the fields are not included or are empty:
%
%       Xsolve,   string     : Pick which Xsolve to use--'FISTA' or
%       'pdNCG.'
%
%       xpos,     bool       :  Constrain X to have nonnegative entries
%           when running XSolve.
%
%       getbias,  bool       : Extract constant bias from observation.
%                                                                                                                      

%% Process input arguments
starttime = tic;
n = size(Y,3);

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

%% PHASE I: First pass at BD
dispfun1 = @(A, X) dispfun(Y, A, X, k, [], 1);

fprintf('PHASE I: \n=========\n');
A = kernel_initialguess;

[A, Xsol, info] = Asolve_Manopt_multikernel( Y, A, lambda1, Xsolve, [], xpos, getbias, dispfun1);
extras.phase1.A = A;
extras.phase1.X = Xsol.X;
extras.phase1.b = Xsol.b;
extras.phase1.info = info;
A_phase_I = A ;
%% PHASE II: Lift the sphere and do lambda continuation
if params.phase2
    k2 = k + 2*kplus;
    dispfun2 = @(A, X) dispfun(Y, A, X, k2, 0, 1);

    A2 = zeros([k2 n]);
    A2(kplus(1)+(1:k(1)), kplus(2)+(1:k(2)), :) = A;
    X2sol = Xsol;
    %X2sol.X = circshift(Xsol.X,-kplus);
    %X2sol.W = circshift(Xsol.W,-kplus);
    % clear A Xsol;

    lambda = lambda1;
    score = zeros(2*kplus+1);
    fprintf('\n\nPHASE II: \n=========\n');
    lam2fac = (lambda2/lambda1)^(1/nrefine);
    i = 1;
    while i <= nrefine + 1
        fprintf('lambda = %.1e: \n', lambda);
        [A2, X2sol, info] = Asolve_Manopt( Y, A2, lambda, Xsolve, X2sol, xpos, getbias, dispfun2 );
        fprintf('\n');

        %Attempt to 'unshift" the a and x by taking the l1-norm over all k-contiguous elements:
        for tau1 = -kplus(1):kplus(1)
            ind1 = tau1+kplus(1)+1;
            for tau2 = -kplus(2):kplus(2)
                ind2 = tau2+kplus(2)+1;
                temp = A2(ind1:(ind1+k(1)-1), ind2:(ind2+k(2))-1,:);
                score(ind1,ind2) = norm(temp(:), 1);
            end
        end
        [temp,ind1] = max(score); [~,ind2] = max(temp);
        tau = [ind1(ind2) ind2]-kplus-1;
        A2 = circshift(A2,-tau);
        X2sol.X = circshift(X2sol.X,tau);
        X2sol.W = circshift(X2sol.W,tau);

        % Save phase 2 extras:
        if i == 1;  idx = 1;    else; idx = i;    end
        extras.phase2(idx).A = A2;
        extras.phase2(idx).X = X2sol.X;
        extras.phase2(idx).b = X2sol.b;
        extras.phase2(idx).info = info;
        if i == 1;  extras.phase2 = fliplr(extras.phase2);  end

        dispfun2(A2,X2sol.X);
        lambda = lambda*lam2fac;
        i = i+1;

    end
end

%% Finished: get the final A, X
if params.phase2
    Aout = A2(kplus(1)+(1:k(1)), kplus(2)+(1:k(2)), :);
    extras.normA = norm(Aout(:));
    Xout = circshift(X2sol.X,kplus) * norm(Aout(:));
    Aout = Aout/norm(Aout(:));
    bout = X2sol.b;
else
    Aout = A;
    extras.normA = norm(Aout(:));
    Xout = Xsol.X;
    bout = Xsol.b;
end

runtime = toc(starttime);
fprintf('\nDone! Runtime = %.2fs. \n\n', runtime);
end

function [ Aout, Xsol, extras ] = Asolve_Manopt_multikernel( Y, Ain, lambda, Xsolve, varargin )
    %ASolve_MANOPT     BD using Manopt solvers.
    %   - Core usage:
    %       [ Aout, Xsol, Stats ] = Asolve_Manopt( Y, Ain, lambda, Xsolve)
    %
    %   - Optional variables:
    %       [ ... ] = Asolve_Manopt( ... , Xinit, Xpos, getbias, dispfun )
    %
    
        load([fileparts(mfilename('fullpath')) '/../config/Asolve_config.mat']); %#ok<*LOAD>
    
        k = size(Ain);
        if (numel(k) > 2)
            num_kernels = k(3); k = k(1:2);
        else
            num_kernels = 1;
        end
        n=1;
        Ain = Ain/norm(Ain(:));
    
        %% Handle the extra variables:
        nvarargin = numel(varargin);
        if nvarargin > 4
            error('Too many input arguments.');
        end
    
        idx = 3;
        if nvarargin < idx || isempty(varargin{idx})
            getbias = false;
        else
            getbias = varargin{idx};
        end
        
        idx = 2;
        if nvarargin < idx || isempty(varargin{idx})
            xpos = false;
        else
            xpos = varargin{idx};
        end
        
        idx = 1;
        if nvarargin < idx || isempty(varargin{idx})
            if strcmp(Xsolve,'FISTA')
                xinit = Xsolve_FISTA(Y, Ain, lambda, mu, [], xpos);
            elseif strcmp(Xsolve,'pdNCG')
                xinit = Xsolve_pdNCG(Y, Ain, lambda, mu, [], xpos);
            end
        else
            xinit = varargin{idx};
        end
    
        idx = 4;
        if nvarargin < idx || isempty(varargin{idx})
            dispfun = @(a, X) 0;
        else
            dispfun = varargin{idx};
        end
    
        %% Set up the problem structure for Manopt and solve
        % The package containing supplement information for the cost, egrad,
        % and ehess functions:
        %{
        suppack.Y = Y;
        suppack.k = k;
        suppack.n = n;
        suppack.mu = mu;
        suppack.xinit = xinit;
        suppack.saveiterates = saveiterates;
        %}
    
        problem.M = spherefactory(prod(k)*n);
        problem.cost = @(a, store) costfun(a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);
        problem.egrad = @(a, store) egradfun(a, store,  Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);
        problem.ehess = @(a, u, store) ehessfun(a, u, store,  Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);
    
        options.statsfun = @(problem, a, stats, store) statsfun( problem, a, stats, store, k, n, saveiterates, dispfun);
        %options.stopfun = @(problem, x, info, last) stopfun(problem, x, info, last, TRTOL);
    
        % Run Manopt solver:
        [Aout, extras.cost, info, extras.options] = trustregions(problem, Ain(:), options);
    
        % Produce final output:
        Aout = reshape(Aout, [k n]);
        if saveiterates
            extras.Aiter = arrayfun(@(i) info(i).A, 1:numel(info), 'UniformOutput', false);
            niter = numel(extras.Aiter);
            extras.Aiter = cell2mat(reshape(extras.Aiter, [1 1 niter]));
            extras.Xiter = arrayfun(@(i) info(i).X, 1:numel(info), 'UniformOutput', false);
            extras.Xiter = cell2mat(reshape(extras.Xiter, [1 1 niter]));
    
            Xsol.X = extras.Xiter(:,:,end);
            Xsol.W = info(end).W;
            Xsol.b = info(end).b;
        else
            if strcmp(Xsolve,'FISTA')
                Xsol = Xsolve_FISTA( Y, Aout, lambda, mu, xinit, xpos, getbias);
            elseif strcmp(Xsolve,'pdNCG')
                Xsol = Xsolve_pdNCG(Y, Aout, lambda, mu, xinit, xpos, getbias);
            end
        end
    end
    
    function [ cost, store ] = costfun( a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve)
        if ~isfield(store, 'X')
            store = computeX( a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);
        end
    
        cost = store.cost;
    end
    
    function [ egrad, store ] = egradfun( a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve)
        if ~isfield(store, 'X')
            store = computeX( a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);
        end
        
        m = size(store.X);
        
        egrad = zeros(prod(k)*n,1);
        
        for i = 1:n
            idx = (i-1)*prod(k) + (1:prod(k));
            tmp = convfft2( store.X, convfft2( reshape(a(idx), k), store.X ) + store.b(i) - Y(:,:,i), 1, m+k-1, m);
            tmp = tmp(1:k(1), 1:k(2));
            egrad(idx) = tmp(:);
        end
    end
    
    function [ ehess, store ] = ehessfun( a, u, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve)
        if ~isfield(store, 'X')
            store = computeX( a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve);
        end
    
        ehess = H_function( u, Y, reshape(a, [k n]), store.X, lambda, mu );
    end
    
    function [ store ] = computeX( a, store, Y, k, n, lambda, mu, xinit, xpos, getbias, Xsolve)
        % Updates the cache to store X*(A), and the active-set whenever a new
        % a new iteration by the trust-region method needs it.
        if strcmp(Xsolve,'FISTA')
            sol = Xsolve_FISTA( Y, reshape(a, [k n]), lambda, mu, xinit, xpos, getbias );
        elseif strcmp(Xsolve,'pdNCG')
            sol = Xsolve_pdNCG( Y, reshape(a, [k n]), lambda, mu, xinit, xpos, getbias);
        end
            
        store.X = sol.X;
        store.W = sol.W;
        store.b = sol.b;
        store.cost = sol.f;
    end
    
    function [ stats ] = statsfun( problem, a, stats, store, k, n, saveiterates, dispfun) %#ok<INUSL>
        if saveiterates
            stats.A = reshape(a, [k n]);
            stats.X = store.X;      % So X could be returned at the end.
            stats.W = store.W;
            stats.b = store.b;
        end
        dispfun(a, store.X);
    end

   function [ Xsol, info ] = Xsolve_FISTA_multikernel( Y, A, lambda, mu, varargin )
    %XSOLVE_FISTA_MULTIKERNEL   Solve for X using FISTA method with multiple kernels
    % Note that this function assumes only one energy layer.
    %   - Core usage:
    %       [ Xsol, info ] = Xsolve_FISTA_multikernel( Y, A, lambda, mu )
    %
    %   - Optional variables:
    %       [ ... ] = Xsolve_FISTA_multikernel( ... , Xinit, Xpos, getbias )
    %       Xinit:      initial value for X
    %       Xpos:       constrain X to be a positive solution
    %       getbias:    extract constant bias as well as X
    %
    
    % Initialize variables and function handles:
    fpath = fileparts(mfilename('fullpath'));
    addpath([fpath '/helpers']);
    load([fpath '/../config/Xsolve_config.mat']); %#ok<*LOAD>
    g = huber(mu);
    num_kernels = size(A,3);
    m = size(Y);
    % define size of X as spatial size + number of kernels
    X_size = [m, num_kernels];

    %% Checking arguments:
    nvararg = numel(varargin);
    if nvararg > 3
        error('Too many input arguments.');
    end

    idx = 1; X = zeros(X_size); b = 0; % note that bias is now a scalar
    if nvararg >= idx && ~isempty(varargin{idx})
        X = varargin{idx}.X;
        b = varargin{idx}.b;
    end

    idx = 2; xpos = false;
    if nvararg >= idx && ~isempty(varargin{idx})
        xpos = varargin{idx};
    end
    
    idx = 3; getbias = false;
    if nvararg >= idx && ~isempty(varargin{idx})
        getbias = varargin{idx};
    end

    % Pre-compute FFT of A for optimization
    A_fft = fft2(A, m(1), m(2));
    R_A = sum(abs(A_fft).^2, 3);
    
    %% Iterate:    
    t=1; W = X; u = b;
    costs = NaN(MAXIT,2);
    doagain = true;  it = 0;  count = 0;
    while doagain
        it = it + 1;
        % Gradients and Hessians:
        grad_fW = zeros(X_size); grad_fu = 0;
        Ri = zeros(m);
        for k = 1:num_kernels
            Ri = Ri + convfft2(A(:,:,k), W(:,:,k));
        end
        Ri = Ri + u - Y;
        grad_fu = sum(Ri(:));

        % Combine loops for optimization
        for k = 1:num_kernels
            grad_fW(:,:,k) = convfft2(A(:,:,k), Ri, 1);
        end

        % FISTA update
        L = max(R_A(:));
        X_ = g.prox(W - 1/L*grad_fW, lambda/L, xpos);
        t_ = (1+sqrt(1+4*t^2))/2;
        W = X_ + (t-1)/t_*(X_-X);
        if getbias
            b_ = u - grad_fu/(2*prod(m));  
            u = b_ + (t-1)/t_*(b_-b);
            b = b_;
        end
        X = X_; t = t_;

        % Check conditions to repeat iteration:
        p = zeros(m);
        for k = 1:num_kernels
            p = p + convfft2(A(:,:,k), X(:,:,k));
        end
        f = norm(p + b - Y, 'fro')^2/2;
        costs(it,1) = f;
        costs(it,2) = g.cost(X, lambda);

        tmp = grad_fW + grad_fu;
        delta = g.diffsubg(X, -tmp, lambda, xpos);
        delta = norm(delta(:))/sqrt(prod(X_size));
        if delta < EPSILON
            count = count+1;
        else
            count = 0;
        end
        doagain = count < 10 && (it < MAXIT);
    end
    
    % Return solution:
    Xsol.X = X;
    Xsol.b = b;
    Xsol.W = W;         % dummy variable for compatibility with pdNCG.
    Xsol.f = sum(costs(it,:));
    info.numit = it;
    info.delta = delta;
    info.costs = costs(1:it,:);
end

function [ out ] = Hxx_function( v, m, A, Hdiag )
    %% 1. Apply rho'*Ca'*P'*P*Ca*rho to v
    % Reshape and extend v
    k = [size(A,1) size(A,2)]; n = size(A,3);
    tmpsz = m + k - 1;
    rhov = extend(reshape(v, m), tmpsz);
    
    % For each slice, apply the rest of the operations
    tmp = zeros(m);
    delta = tmpsz - m;
    mask = false(tmpsz); 
    mask(floor(delta(1)/2)+(1:m(1)), floor(delta(2)/2)+(1:m(2))) = true;
    for i = 1:n
        % apply P'*P*Ca to rhov
        currslice = mask.*cconvfft2(A(:,:,i),rhov);
        
        % apply rho'*Ca' and save
        currslice = extend(cconvfft2(A(:,:,i),currslice,tmpsz,'left'),m,1);
        tmp = tmp + currslice;
    end
    out = tmp(:) + Hdiag.*v;
end

function [ h ] = huber( mu )
    if nargin < 1 || isempty(mu)
        mu = 1e-6;
    end

    h.cost = @(x, lambda) cost(x, lambda, mu);
    h.prox = @(x, lambda, xpos) prox (x, lambda, mu, xpos);
    h.diffsubg = @(x, y, lambda, xpos) diffsubg(x, y, lambda, mu, xpos);
end

function [ hx ] = cost(x, lambda, mu)
    leq = abs(x) <= mu;
    hx = sum(x(leq).^2/2/mu) + sum(abs(x(~leq)) - mu/2);
    hx = lambda*hx;
end

function [ proxh ] = prox(x, lambda, mu, xpos)
    if nargin < 4 || isempty(xpos)
        xpos = false;
    end

    leq = abs(x) <= lambda + mu;
    proxh = zeros(size(x));
    proxh(leq) = x(leq)./(1+lambda/mu);
    proxh(~leq) = x(~leq) - lambda*sign(x(~leq));

    if xpos;  proxh = max(proxh,0);  end
end

function [ eps ] = diffsubg(x, y, lambda, mu, xpos)
%DIFFSUBG  Difference from y to the subgradient of the loss function at x.
%  In this case the Huber loss function is smooth so only need to subtract
%  from the gradient.

    leq = abs(x) <= mu;
    subg = zeros(size(x));

    subg(leq) = x(leq)/mu;
    subg(~leq) = sign(x(~leq));
    if xpos;  subg(x<0) = 0;  end

    eps = y - lambda*subg;
end

function [ H_v ] = H_function( v, Y, A, X, lambda, mu )
%H_FUNCTION     Apply the Euclidean Hessian to a vector.
    load([fileparts(mfilename('fullpath')) '/../../config/Hfunction_config.mat']);
    
    m = size(Y); k = size(A); n = size(Y,3);
    m = m(1:2); k = k(1:2);
    tmpsz = m + k - 1;
    
    Haa_v = zeros([k n]);
    r = zeros([m n]);
    Hxa_v = zeros(m);
    for i = 1:n
        idx = (i-1)*prod(k) + (1:prod(k));
        vi = reshape(v(idx), k);
        Haa_v(:,:,i) = convfft2( X, convfft2(X, vi, 0, tmpsz, m), 1, tmpsz, k);
    
        r(:,:,i) = convfft2(A(:,:,i),X) - Y(:,:,i);
        Hxa_v = Hxa_v ...
            + convfft2( A(:,:,i), convfft2(X, vi, 0, tmpsz, m), 1) ...
            + Hxres( r(:,:,i), vi, k, m );
    end
    
    hesspendiag = lambda * mu^2*(mu^2 + X(:).^2).^(-3/2);    
    Hxxfun = @(u) Hxx_function(u, m, A, hesspendiag);
    pcgprecond = @(u) u./(1 + hesspendiag);
    [HxxInv_Hxa_v,~] = pcg(Hxxfun, Hxa_v(:), PCGTOL, PCGIT, pcgprecond);
    HxxInv_Hxa_v = reshape(HxxInv_Hxa_v, m);
    
    % Hax = Ikm'*(CA*CX' + (CA*CX-CY)*Pi);
    
    Hax_HxxInv_Hxa_v = zeros([k n]);
    for i = 1:n
        Hax_HxxInv_Hxa_v(:,:,i) = ...
            convfft2( X, convfft2(A(:,:,i), HxxInv_Hxa_v, 0, tmpsz, m), 1, tmpsz, k) ...
            + Hxres( r(:,:,i), HxxInv_Hxa_v, m, k );
    end
    
    H_v = Haa_v - Hax_HxxInv_Hxa_v;
    H_v = H_v(:);

end

function [ out ] = Hxres( res, in, insz, outsz )
    tmpsz = outsz + insz - 1;
    Vhat = fft2(extend(reshape(in, insz), tmpsz));
    reshat = fft2(shrink(res, tmpsz, 1));
    out = extend(ifft2(reshat.*conj(Vhat)), outsz, 1 );
end