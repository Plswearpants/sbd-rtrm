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

%% Phase I: Initialization and First Iteration
fprintf('PHASE I: Initialization and First Iteration\n');

A = kernel_initialguess;
X_struct = struct();
Xiter = zeros([size(Y),kernel_num]);
biter = zeros(kernel_num,1);
kernel_quality_factors = zeros(maxIT,kernel_num);
observation_quality_factors = zeros(maxIT, 1);

for iter = 1:maxIT
    starttime = tic;

    % Compute Y_background
    Y_background = Y;
    for m = 1:kernel_num
        if iter > 1
            Y_background = Y_background - convfft2(A{m}, Xiter(:,:,m));
        end
    end

    % Update each kernel
    parfor n = 1:kernel_num
        Yiter = Y_background + (iter > 1) * convfft2(A{n}, Xiter(:,:,n));
        dispfun1 = @(A, X) dispfun{n}(Y, A, X, k(n,:), []);
        
        if iter == 1
            % Initial X computation (previously in Phase 0)
            X_struct.(['x',num2str(n)]) = Xsolve_FISTA_tunable(Y, A{n}, lambda1(n), mu, AX_iteration, [], xpos);
        end
        
        [A{n}, X_struct.(['x',num2str(n)]), info] = Asolve_Manopt_tunable(Yiter, A{n}, lambda1(n), Xsolve, AX_iteration, X_struct.(['x',num2str(n)]), xpos, getbias, dispfun1);
        
        Xiter(:,:,n) = X_struct.(['x',num2str(n)]).X;
        biter(n) = X_struct.(['x',num2str(n)]).b;
        
        % Compute kernel quality factor
        kernel_quality_factors(iter, n) = compute_kernel_quality{n}(A{n});
    end

    % Compute observation quality factor
    Yiter_combined = sum(Yiter, 3) - (kernel_num-1)*Y_background;
    observation_quality_factors(iter) = var(Y(:) - Yiter_combined(:));

    runtime = toc(starttime);
    fprintf('Iteration %d: Runtime = %.2fs, Observation Quality Factor = %.8e\n', iter, runtime, observation_quality_factors(iter));
    for n = 1:kernel_num
        fprintf('Kernel %d Quality Factor = %.8e\n', n, kernel_quality_factors(iter, n));
    end
end

% Store results and quality factors
extras.observation_quality_factors = observation_quality_factors;
extras.kernel_quality_factors = kernel_quality_factors;
extras.phase1.Aout = A;
extras.phase1.Xout = Xiter;
extras.phase1.biter = biter;

% Visualize quality factors
figure;
subplot(2,1,1);
semilogy(1:maxIT, observation_quality_factors);
xlabel('Iteration');
ylabel('Observation Quality Factor');
title('Observation Quality Factor vs. Iteration');
grid on;

subplot(2,1,2);
hold on;
for n = 1:kernel_num
    semilogy(1:maxIT, kernel_quality_factors(:,n), 'DisplayName', sprintf('Kernel %d', n));
end
xlabel('Iteration');
ylabel('Kernel Quality Factor');
title('Kernel Quality Factors vs. Iteration');
legend('show');
grid on;
hold off;

%% PHASE II: Lift the sphere and do lambda continuation
if params.phase2
    fprintf('\n\nPHASE II: \n=========\n');
    k2 = k + 2*kplus;

    A2 = cell(1, kernel_num);
    X2_struct = struct();
    for n = 1:kernel_num
        A2{n} = zeros(k3(n,:));
        A2{n}(kplus(n,1)+(1:k(n,1)), kplus(n,2)+(1:k(n,2))) = A{n};
        X2_struct.(['x',num2str(n)]) = X_struct.(['x',num2str(n)]);
        % X2_struct.(['x',num2str(n)]).X = circshift(X_struct.(['x',num2str(n)]).X, -kplus(n,:));
        % X2_struct.(['x',num2str(n)]).W = circshift(X_struct.(['x',num2str(n)]).W, -kplus(n,:));
    end

    lambda = lambda1;
   
    lam2fac = (lambda2./lambda1).^(1/nrefine);
    
    for i = 1:nrefine + 1
        fprintf('lambda iteration %d/%d: \n', i, nrefine + 1);
        
        % initial Y_background
        Y_background = Y;
        for m = 1:kernel_num
            Y_background = Y_background - convfft2(A2{m}, X2_struct.(['x',num2str(m)]).X);
        end

        for n = 1:kernel_num
            fprintf('Processing kernel %d, lambda = %.1e: \n', n, lambda(n));
            Yiter = Y_background + convfft2(A2{n}, X2_struct.(['x',num2str(n)]).X);
            dispfun2 = @(A, X) dispfun{n}(Y, A, X, k2(n,:), 0, 1);
            [A2{n}, X2_struct.(['x',num2str(n)]), info] = Asolve_Manopt_tunable(Yiter, A2{n}, lambda(n), Xsolve, AX_iteration, X2_struct.(['x',num2str(n)]), xpos, getbias, dispfun2);
            
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

            % Save phase 2 extras:
            extras.phase2.A{n} = A2{n};
            extras.phase2.X{n} = X2_struct.(['x',num2str(n)]).X;
            extras.phase2.b(n) = X2_struct.(['x',num2str(n)]).b;
            extras.phase2.info{n} = info;
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

end
