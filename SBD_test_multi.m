function [ Aout, Xout, bout, extras ] = SBD_test_multi( Y, k, params, dispfun, kernel_initialguess, AX_iteration, maxIT)
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
compute_kernel_quality = params.compute_kernel_quality;  % Extract the function handle

%% PHASE I: 1st iteration through BD

fprintf('PHASE 0: \n=========\n');
A = kernel_initialguess;

Xint = zeros([size(Y),kernel_num]); % preallocate memory for Xint
Xiter = zeros([size(Y),kernel_num]); % preallocate memory for Xiter, the intermediate X that will be updated iteratively
Yiter = zeros([size(Y),kernel_num]); % preallocate memory for Yiter, the intermediate Y that will be updated iteratively
biter = zeros(kernel_num,1); % preallocate memory for biter, the bias term that will be updated iteratively
Kernel_quality_factors = zeros(maxIT,kernel_num);
Observation_quality_factors = zeros(maxIT, 1);

% initial Xint for each kernel
X_struct = struct();

parfor n = 1 : kernel_num
    X_struct.(['x',num2str(n)])  = Xsolve_FISTA_tunable(Y, A{n}, lambda1(n), mu, AX_iteration,[], xpos);
    Xint(:,:,n) = X_struct.(['x',num2str(n)]).X;   
end
fprintf('finished initializing Xint\n');

% initial Y_background, which is the residual after removing the contribution of each kernel
for m = 1: kernel_num
    % update Y_background
    Y_background = (Y - convfft2(A{m}, Xint(:,:,m)));
end

% iteratively update each kernel
parfor n = 1:kernel_num
    Yiter(:,:,n) = Y_background + convfft2(A{n}, Xint(:,:,n));
    dispfun1 = @(A, X) dispfun{n}(Y, A, X, k(n,:), []);
    [A{n}, X_struct.(['x',num2str(n)]), info] = Asolve_Manopt_tunable( Yiter(:,:,n), A{n}, lambda1(n), Xsolve, AX_iteration, X_struct.(['x',num2str(n)]), xpos, getbias, dispfun1);
    Xiter(:,:,n) = X_struct.(['x',num2str(n)]).X;
    % update Yiter
    Yiter(:,:,n) = Y_background + convfft2(A{n}, Xiter(:,:,n));
end 
fprintf('finished initializing A\n');

Yiter_initial = sum(Yiter,3)-(kernel_num-1)*Y_background; % Yiter_initial resembles the initial Y after this iteration, with the excessive background removed
Observation_quality_factor = var(Y(:)-Yiter_initial(:)); % the quality factor is the variance of the residual

% Compute initial kernel quality factors
for n = 1:kernel_num
    Kernel_quality_factors(1, n) = compute_kernel_quality{n}(A{n});
end

%% Phase II : Use the first iteration result and Pass to BD for batch iteration
fprintf('PHASE II: \n=========\n');

Observation_quality_factors(1) = Observation_quality_factor;  % Store the initial quality factor

for iter = 2:maxIT
    starttime = tic;

    % initial Y_background, which is the residual after removing the contribution of each kernel
    Y_background = Y;
    for m = 1:kernel_num
        Y_background = Y_background - convfft2(A{m}, Xiter(:,:,m));
    end

    % iteratively update each kernel
    parfor n = 1:kernel_num
        Yiter(:,:,n) = Y_background + convfft2(A{n}, Xiter(:,:,n));
        dispfun1 = @(A, X) dispfun{n}(Y, A, X, k(n,:), []);
        [A{n}, X_struct.(['x',num2str(n)]), info] = Asolve_Manopt_tunable(Yiter(:,:,n), A{n}, lambda1(n), Xsolve, AX_iteration, X_struct.(['x',num2str(n)]), xpos, getbias, dispfun1);
        % update Xiter
        Xiter(:,:,n) = X_struct.(['x',num2str(n)]).X;
        biter(n) = X_struct.(['x',num2str(n)]).b;
        % update Yiter
        Yiter(:,:,n) = Y_background + convfft2(A{n}, Xiter(:,:,n));
        
        % Compute kernel quality factor
        Kernel_quality_factors(iter, n) = compute_kernel_quality{n}(A{n});
    end 

    Yiter_this = sum(Yiter, 3)-(kernel_num-1)*Y_background;
    Observation_quality_factors(iter) = var(Y(:) - Yiter_this(:));
    runtime = toc(starttime);
    fprintf('Iteration %d: Runtime = %.2fs, Observation Quality Factor = %.8e\n', iter, runtime, Observation_quality_factors(iter));
    for n = 1:kernel_num
        fprintf('Kernel %d Quality Factor = %.8e\n', n, Kernel_quality_factors(iter, n));
    end
end

% Store quality factors in extras
extras.Observation_quality_factors = Observation_quality_factors;
extras.Kernel_quality_factors = Kernel_quality_factors;

%% Finished: get the final A, X
Aout = A;
Xout = Xiter;
bout = biter;

% Visualize quality factors
figure;
subplot(2,1,1);
semilogy(1:maxIT, Observation_quality_factors);
xlabel('Iteration');
ylabel('Observation Quality Factor');
title('Observation Quality Factor vs. Iteration');
grid on;

subplot(2,1,2);
hold on;
for n = 1:kernel_num
    semilogy(1:maxIT, Kernel_quality_factors(:,n), 'DisplayName', sprintf('Kernel %d', n));
end
xlabel('Iteration');
ylabel('Kernel Quality Factor');
title('Kernel Quality Factors vs. Iteration');
legend('show');
grid on;
hold off;

%{
extras.phase1.A = A;
extras.phase1.X = Xsol.X;
extras.phase1.b = Xsol.b;
extras.phase1.info = info;
A_phase_I = A ;


%% PHASE III: Lift the sphere and do lambda continuation
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
%}
end
