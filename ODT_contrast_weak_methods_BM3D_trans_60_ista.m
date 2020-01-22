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

scatteringStrength = 1;

bm3dTausFull = [0.001494697762586, 0.002075346882639, 0.001943857936476, ...
                0.002065028404822, 0.001825286916654, 0.001915072490250, ...
                0.001833798964586, 0.001873208929149];
            
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
stats = cell(num, length(batchSize));
for i = 1:num
    %% generate data
    truth = truths(:,:,i);
    f = truth*scatteringStrength;
 
    % signal-to-noise ration in dB
    inputSNR = 40;

    %%% measurements
%     utotDom = forwardProp(em.uincDom, f, em.domainGreensFunction, em.uincDom,...
%         em.dx, em.dy);  use linear scattering
    uscatPred = fullPropagateToSensor(f, em.uincDom,...
        em.sensorGreensFunction, em.receiverMask, em.dx, em.dy);
    data = addawgn(uscatPred, inputSNR);

    for j = 1:length(batchSize)
        %% common reconstruction parameters
        accelerate = false;
        verbose = true;
        plotRecon = false;
        amutTime = inf;
        numIter = 2e3;
        step = 3e2;
        stochastic = true;
        keepNum = batchSize(j);

        % tv parameters
        tvProxIter = 20;
        bounds = [0, inf];

        %% first-Born and BM3D

        % optimize only for the first image
        tau = bm3dTausFull(i);

        dobj = FirstBornClass(data, em);
        robj = BM3DClass(imgSize, tau);

        %%% Stochastic first-Born and BM3D
        fprintf('reconstruct Stochastic first-Born and BM3D\n');

        % reconstruct
        [recon, outs] = fistaEst(dobj, robj,...
            'accelerate', accelerate,...
            'verbose', verbose,...
            'amuttime',amutTime,...
            'plotRecon', plotRecon,...
            'xtrue', f,...
            'numIter', numIter,...
            'step', step,...
            'stochastic', stochastic,...
            'keepnum', keepNum);

        % store recons and stats
        outs.recon = recon;
        outs.tau   = tau;
        outs.keepNum = keepNum;
        stats{i,j} = outs;
    end
end

if(accelerate)
    save('ODT_fista_contrast=weak_methods=bm3d_trans=60.mat', 'stats');
else
    save('ODT_ista_contrast=weak_methods=bm3d_trans=60.mat', 'stats');
end
