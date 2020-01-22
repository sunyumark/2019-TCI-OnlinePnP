function [uhat, outs] = computeTotalField(uin, f, g, dx, dy, varargin)
%%% Computes the total field inside the object given the input field,
%%% scattering potential, and the Green's function.
%%%
%%% Input:
%%% - uin: (Ny x Nx) input field
%%% - f: (Ny x Nx) scattering potential
%%% - g: (2*Ny x 2*Nx) Green's function
%%% - dx, dy: sampling steps
%%%
%%% Ouput:
%%% - uhat: (Ny x Nx) total field
%%%
%%% U. S. Kamilov, MERL, 2016.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ny, Nx] = size(uin);

numIter = 100;
plotRecon = false;
step = 1e-3;
tol = 1e-4;
verbose = false;
uhat0 = uin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nargs = length(varargin); % Number of options

for i = 1:2:nargs % Go through options
    
    name = lower(varargin{i}); % Extract name
    value = varargin{i+1}; % Extract value
    
    switch(name)
        case 'numiter'
            numIter = value;
        case 'plotrecon'
            plotRecon = value;
        case 'step'
            step = value;
        case 'tol'
            tol = value;
        case 'verbose'
            verbose = value;
        case 'uhat0'
            uhat0 = value;
        otherwise
            error('computeTotalField: input is not recognized!');
    end
end

assert(all((2*size(f))==size(g)), 'all((2*size(f))==size(g))');
assert(all(size(f)==[Ny, Nx]), 'all(size(f)==[Ny, Nx])');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% useful handles
A = @(u) A_forward(u, f, g, dx, dy);
AT = @(z) A_adjoint(z, f, g, dx, dy);
evaluateCost = @(u) sum(sum(abs(uin-A(u)).^2))/sum(sum(abs(uin).^2));
evaluateTol = @(u,unext) norm(u(:)-unext(:))/norm(u(:));

%%% Initialize
uhat = uhat0;
s = uhat;
q = 1;

%%% Track evolution
outs.time = zeros(numIter, 1);
outs.cost = zeros(numIter, 1);
outs.relativeChange = zeros(numIter, 1);

if(plotRecon)
    figureHandle = figure('Name', 'computeTotalField');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for indIter = 1:numIter
    timeStart = cputime;
    
    uhatnext = s-step*AT(A(s)-uin);
    
    qnext = 0.5*(1+sqrt(1+4*q*q));
    
    s = uhatnext + ((q-1)/qnext)*(uhatnext-uhat);
    
    outs.time(indIter) = cputime - timeStart;
    outs.cost(indIter) = evaluateCost(uhatnext);
    outs.relativeChange(indIter) = evaluateTol(uhat, uhatnext);
    
    q = qnext;
    uhat = uhatnext;
    
    if(verbose)
        fprintf('[computeTotalField: %d/%d]', indIter, numIter);
        fprintf('[cost: %.2e]', outs.cost(indIter));
        fprintf('[time: %.1f]', sum(outs.time));
        fprintf('[tols: %.1e]', outs.relativeChange(indIter));
        fprintf('\n');
    end
    
    if(plotRecon)
        figure(figureHandle);
        set(figureHandle, 'Color', 'w',...
            'Name',sprintf('computeTotalField [%d/%d]',indIter,numIter));
        
        subplot(2, 2, 1);
        imagesc(abs(uhat(:,:)./uin(:,:)));
        axis equal off;
        title('abs-uscat');
        set(gca, 'FontSize', 16);
        colorbar;
        
        subplot(2, 2, 2);
        imagesc(angle(uhat(:,:)./uin(:,:)));
        axis equal off;
        title('angle-uscat');
        set(gca, 'FontSize', 16);
        colorbar;
        
        subplot(2, 2, 3);
        semilogy(1:indIter, outs.cost(1:indIter), 'b-',...
            indIter, outs.cost(indIter), 'ro', 'LineWidth', 1.5);
        xlim([1 numIter]);
        grid on;
        title(sprintf('cost: %.2e', outs.cost(indIter)));
        set(gca, 'FontSize', 16);
        
        subplot(2, 2, 4);
        semilogy(1:indIter, outs.relativeChange(1:indIter), 'b-',...
            indIter, outs.relativeChange(indIter), 'ro', 'LineWidth', 1.5);
        xlim([1 numIter]);
        grid on;
        title(sprintf('relativeChange: %.1e', outs.relativeChange(indIter)));
        set(gca, 'FontSize', 16);
        
        drawnow;
    end
    
    if(outs.relativeChange(indIter) < tol)
        outs.cost = outs.cost(1:indIter);
        outs.relativeChange= outs.relativeChange(1:indIter);
        outs.time = outs.time(1:indIter);
        break;
    end
end