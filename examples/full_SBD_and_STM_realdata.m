clc; clear;

% Import Manopt and initialize the SBD package
run('../init_sbd');
fprintf('\n\n');

%% 0. Load the .3ds data

% INPUTS
% 1: Data file to load, including file type ('QPI.3ds' for example)
% 2: Smoothing sigma for current data

% OUTPUTS
% header: Variable containing all experimental parameters
% I: Current data, smoothed by sigma
% dIdV: Numerically differentiated current data
% voltage: Vector of voltages for current
% midV: Vector on voltages for dIdV/QPI (midpoint of voltage vector)
% QPI: Fourier transformed dIdV data

% Modified function load3dsall from supplied matlab code from Nanonis
[header, par, I, dIdV, LockindIdV, bias, midV, QPI, LockinQPI] = load3dsall('Grid_20210906_HR_CZ.3ds', 5);
xsize = header.grid_dim(1);
ysize = header.grid_dim(2);
elayer = header.points;
estart = par(1);
eend = par(2);
%% 0.1 Slice selection
d3gridDisplay(dIdV,'dynamic');
selected_slice = input('Enter the slice number you want to analyze: ');
%% I. SIMULATE DATA FOR SBD:
%  =========================

%% 1. Kernel settings - see kernel types below
kerneltype = 'simulated_STM';   
n = 1;               	% number of kernel slices
k = squareDrawSize(dIdV(:,:,selected_slice));           	% determine kernel size

%% 2. Activation map generation:
% Generate activation map based on the sliced data
X0=activationCreateClick(dIdV(:,:,40));

m = size(X0);          % image size for each slice / observation grid

eta = estimate_noise(dIdV(:,:,selected_slice));             % additive noise variance

%% 3. Generate kernel
switch kerneltype
    case 'random'
    % Randomly generate n kernel slices
        A0 = randn([k n]);
    
    case 'simulated_STM'
    % Randomly choose n kernel slices from simulated LDoS data
        load('example_data/LDoS_sim.mat');
        sliceidx = randperm(size(LDoS_sim,3), n);
        
        A0 = NaN([k n]);
        for i = 1:n
            A0 = imresize(LDoS_sim(:,:,sliceidx), k);
        end
        
    otherwise
        error('Invalid kernel type specified.')
end

% Need to put each slice back onto the sphere
A0 = proj2oblique(A0);

%% 5 observation generation:
Y = zeros([m n]);
for i = 1:n                           	% observation
    Y(:,:,i) = convfft2(A0(:,:,i), X0);     
end
Y = Y + sqrt(eta)*randn([m n]);

%% II. Sparse Blind Deconvolution:
%  ===============================
%% 1. Settings - refer to documents on details for setting parameters.

% A function for showing updates as RTRM runs
dispfun = @( Y, A, X, k, kplus, idx ) showims(Y,A0,X0,A,X,k,kplus,idx);

% SBD settings
params.lambda1 = 1e-1;              % regularization parameter for Phase I

params.phase2 = true;               % whether to do Phase II (refinement)
params.kplus = ceil(0.2 * k);       % padding for sphere lifting
params.lambda2 = 5e-2;              % FINAL reg. param. value for Phase II
params.nrefine = 3;                 % number of refinements

% Want entries of X to be nonnegative: see SBD_main.m
params.signflip = 0.2;
params.xpos     = true;
params.getbias  = true;
params.Xsolve = 'FISTA';

% 2. The fun part
[Aout, Xout, extras] = SBD( Y, k, params, dispfun );
