function [x, outs] = fistaEst(dObj, rObj, varargin)
% Generic object oriented implementation of FISTA
%
% Input:
% - dObj: data-fidelity term
% - rObj: regularization term
%   - varargin: options for the algorithm
%     Options:
%     - backTol: tolerance of backtracking (def: 1e-12)
%     - eta: backtracking constant (def: 1/2)
%     - numIter: number of iterations (def: 100)
%     - plotRecon: plot reconstruction (def: false)
%     - step: step-size for FISTA (def: 1)
%     - tol: stopping criterion for the algorithm (def: 1e-4)
%     - tolCount: number of satisfactions of tol before stop (def: 3)
%     - verbose: print command line message (def: false)
%     - xtrue: oracle (i.e. the original signal)
%     - xinit: initialization
%     - stochastic: if use stochastic version or not false;
%     - keepNum: num of measurements the alg use for each iteration;
%
% Output:
% - x: estimated image
% - outs: extra data
%     - outs.obj: per iteration obj
%     - outs.snr: per iteration snr
%     - outs.time: per iteration CPU-time
%     - outs.relativeChange: per iteration distance from tolerance
%
% U. S. Kamilov, CIG, WUSTL, 2017.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigSize = dObj.size;

accelerate = true;
numIter = 100;
amutTime = inf;
plotRecon = false;
step = 1;
tol = 1e-20;
tolCount = 3;
verbose = false;
xinit = zeros(sigSize);
stochastic = false;
keepNum = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nargs = length(varargin); % Number of options

for i = 1:2:nargs % Go through options
    
    name = lower(varargin{i}); % Extract name
    value = varargin{i+1}; % Extract value
    
    switch(name)
        case 'accelerate'
            accelerate = value;
        case 'numiter'
            numIter = value;
        case 'plotrecon'
            plotRecon = value;
        case 'step'
            step = value;
        case 'tol'
            tol = value;
        case 'tolcount'
            tolCount = value;
        case 'verbose'
            verbose = value;
        case 'xtrue'
            xtrue = value;
        case 'xinit'
            xinit = value;
        case 'stochastic'
            stochastic = value;
        case 'keepnum'
            keepNum = value;
        case 'amuttime'
            amutTime = value;
        otherwise
            error('fistaEst: input is not recognized!');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xtrueSet = exist('xtrue', 'var'); % true signal is set

% useful handles
evaluateSnr = @(xtrue, x) 20*log10(norm(xtrue(:))/norm(xtrue(:)-x(:)));
evaluateTol = @(x,xnext) norm(x(:)-xnext(:))/norm(x(:));

% initialize iterates
x = xinit;
s = x;
t = 1;

% initialize measurements mask
totNum = size(dObj.y, 1); % extract number of measurements
id = 0;
keepIdx = zeros(1,keepNum);
for j = 1:keepNum
    keepIdx(j) = floor(id);
    id = id + totNum/keepNum;
end

% initializer for proximal
p = rObj.init;
pfull = rObj.init;

% diagnostic data
outs.time = zeros(numIter, 1);
outs.obj = zeros(numIter, 1);
outs.relativeChange = zeros(numIter, 1);
if(xtrueSet)
    outs.snr = zeros(numIter, 1);
end

if(plotRecon)
    h = figure('Name', 'fistaEst');
end

tolCounter = tolCount; % conter for tolerance violations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for indIter = 1:numIter
    
    timeStart = cputime;    
    % get: d = D(s) and g = grad(D(s))
    % tentative update (pnext = trial dual variable)
    if(stochastic)
        [g, ~, keepIdx] = dObj.gradSub(s, keepIdx);
        [xnext, pnext] = rObj.prox(s-step*g, step, p, 1);
    else
        g = dObj.grad(s);
        [xnext, pnext] = rObj.prox(s-step*g, step, p, 1);
    end
    timeStop = cputime;
    
    % calculate full gradient for objective function
    gfull = dObj.grad(x);
    [Px, pfull] = rObj.prox(x-step*gfull, step, pfull, 1);
    
    % store obj0
    if indIter == 1
       outs.obj0 = norm(x(:)-Px(:))^2;  % should be vector norm
    end
    
    % store objective value
    outs.obj(indIter) = norm(x(:)-Px(:))^2/outs.obj0; % should be vector norm
    
    % acceleration
    if(accelerate)
        tnext = 0.5*(1+sqrt(1+4*t*t));
    else
        tnext = 1;
    end
    s = xnext + ((t-1)/tnext)*(xnext-x);

    % diagnostic data
    outs.time(indIter) = timeStop - timeStart;
    outs.relativeChange(indIter) = evaluateTol(x, xnext);
    if(xtrueSet)
        outs.snr(indIter) = evaluateSnr(xtrue, x);
    end
    
    % update iterates
    p = pnext;
    t = tnext;
    x = xnext;  
    
    % print to command line
    if(verbose)
        fprintf('[fistaEst: %d/%d]', indIter, numIter);
        fprintf('[tols: %.5e]', outs.relativeChange(indIter));
        fprintf('[step: %.1e]', step);
        if(xtrueSet)
            fprintf('[snr: %.2f]', outs.snr(indIter));
        end
        fprintf('[||x-P(x)||^2: %.3e]', outs.obj(indIter));
        fprintf('[time: %.1f]', sum(outs.time));
        fprintf('\n');
    end
    
    % plot graphically
    if(plotRecon)
        figure(h);
        set(h, 'Name', sprintf('fistaEst [%d/%d]', indIter, numIter));
        
        subplot(3, 2, 1:4);
        dObj.draw(x);
        
        if(~xtrueSet)
            subplot(3, 2, 5:6);
            semilogy(1:indIter, outs.obj(1:indIter), 'b-',...
                indIter, outs.obj(indIter), 'ro', 'LineWidth', 1.5);
            xlim([1 numIter]);
            grid on;
            title(sprintf('||x-P(x)||^2: %.4e', outs.obj(indIter)));
            set(gca, 'FontSize', 18);
        else
            subplot(3, 2, 5);
            semilogy(1:indIter, outs.obj(1:indIter), 'b-',...
                indIter, outs.obj(indIter), 'ro', 'LineWidth', 1.5);
            xlim([1 numIter]);
            grid on;
            title(sprintf('||x-P(x)||^2: %.4f', outs.obj(indIter)));
            set(gca, 'FontSize', 18);
            
            subplot(3, 2, 6);
            plot(1:indIter, outs.snr(1:indIter), 'b-',...
                indIter, outs.snr(indIter), 'ro', 'LineWidth', 1.5);
            xlim([1 numIter]);
            grid on;
            title(sprintf('snr: %.2f', outs.snr(indIter)));
            set(gca, 'FontSize', 18);
        end
        drawnow;
    end
    
    % stopping criterion
    if(tolCounter <= 0)
        if(xtrueSet)
            outs.snr = outs.snr(1:indIter);
        end
        outs.obj = outs.obj(1:indIter);
        outs.relativeChange= outs.relativeChange(1:indIter);
        outs.time = outs.time(1:indIter);
        break;
    end
    
    if(sum(outs.time(:)) > amutTime)
        if(xtrueSet)
            outs.snr = outs.snr(1:indIter);
            outs.dist = outs.dist(1:indIter);
        end
        outs.obj  = outs.obj(1:indIter);
        outs.relativeChange= outs.relativeChange(1:indIter);
        outs.time = outs.time(1:indIter);
        break;
    end
        
    % check relative change
    if(outs.relativeChange(indIter) < tol)
        tolCounter = tolCounter - 1;
    else
        tolCounter = tolCount;
    end
    
end