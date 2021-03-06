%%% This script generates Cisor reconstruction
%%%
%%% Yu Sun, CIG, WUSTL, 2017.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; home;

addpath(genpath(pwd));

rng('default'); % test seeds is equal to 128

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scattering contrast [strength of scattering]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scatteringStrength = 1; % 1E-6; % try 1e-1 and 1e-6

tvTausFull = [3.045837411983810e-07, 5.529754415521786e-07, 9.728730768317722e-07, ...
              7.288032332562928e-07, 4.953506599971617e-07, 4.882124587925646e-07, ...
              6.074982097753717e-07, 5.280238095238872e-07];
          
batchSize = [10,20,30];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loaded = load('testImages');
truths = loaded.X;

% size of the truths
[nx,ny,num] = size(truths);
imgSize = [nx,ny];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% shortcuts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

evaluateSNR = @(x, xhat) 20*log10(norm(x(:))/norm(x(:)-xhat(:)));
evaluatePSNR = @(x, xhat) 20*log10(1) - 10*log10(norm(x(:)-xhat(:))^2 / (nx*ny));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% em and discretization parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% size of the reconstruction domain
em.Lx = 0.18; % [m]
em.Ly = 0.18; % [m]

% number of pixels
em.Nx = imgSize(1);
em.Ny = imgSize(2);

% smallest distance between objects
em.dx = em.Lx/em.Nx; % [m]
em.dy = em.Lx/em.Ny; % [m]

% locations of the pixels
em.x = (-em.Nx/2:em.Nx/2-1)*em.dx;
em.y = (-em.Ny/2:em.Ny/2-1)*em.dy;

% speed of light [m/s]
em.c = 299792458;

% wavelength [m]
em.lambda = em.dx*6;

% measured frequency [GHz]
em.freq = em.c./em.lambda/1e9;

% number of receivers
em.numRec = 360;

% number of transmissions
em.numTrans = 60;

% radius of a rig where sensors are located [m]
em.sensorRadius = 1.6;

% wavenumber [1/m]
em.kb = 2*pi./em.lambda;

% generate useful functions
em = generateEmFunctions(em);

%% the for loop
stats = cell(num,length(batchSize));
for i = 1:num
    %% generate data
    truth = truths(:,:,i);
    f = truth*scatteringStrength;
 
    % signal-to-noise ration in dB
    inputSNR = 40;

    %%% measurements
    uscatPred = fullPropagateToSensor(f, em.uincDom,...
        em.sensorGreensFunction, em.receiverMask, em.dx, em.dy);
    data = addawgn(uscatPred, inputSNR);

    for j = 1:length(batchSize)
    
        %% common reconstruction parameters
        accelerate = false;
        verbose = true;
        plotRecon = true;
        amutTime = inf;
        numIter = 2e3;
        step = 3e2; % The largest step size by using backtracking
        stochastic = true;
        keepNum = batchSize(j);

        % tv parameters
        tvProxIter = 20;
        bounds = [0, inf];

        %% first-Born and total variation

        %%% reconstruction
        fprintf('reconstruct Stochastic first-Born and total variation\n');

        % optimize only for the first image
        tau = tvTausFull(i);

        dobj = FirstBornClass(data, em);
        robj = TVClass(imgSize, tau, 'maxiter', tvProxIter, 'bounds', bounds);

        % reconstruct
        [recon, outs] = fistaEst(dobj, robj,...
            'accelerate', accelerate,...
            'verbose', verbose,...
            'plotRecon', plotRecon,...
            'xtrue', f,...
            'numIter', numIter,...
            'amuttime', amutTime,...
            'step', step,...
            'stochastic', stochastic,...
            'keepnum', keepNum);

        % store recons and stats
        outs.recon = recon;
        outs.tau   = tau;
        outs.keepNum = keepNum;
	outs.step = step;
        stats{i,j} = outs;
    end
end

if(accelerate)
    save('ODT_fista_contrast=weak_methods=tv_trans=60.mat', 'stats');
else
    save('ODT_ista_contrast=weak_methods=tv_trans=60.mat', 'stats');
end
