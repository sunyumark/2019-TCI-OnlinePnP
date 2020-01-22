function [utotDom, outs] = pcgWrap(A, b, varargin)
%%% Wrapper function for MATLAB's pcg
%%%
%%% U. S. Kamilov, MERL, 2017


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ny, Nx] = size(b);

numIter = 1000;
plotRecon = false;
tol = 1e-6;
xinit = zeros(Ny, Nx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nargs = length(varargin); % Number of options

for i = 1:2:nargs % Go through options
    
    name = lower(varargin{i}); % Extract name
    value = varargin{i+1}; % Extract value
    
    switch(name)
        case 'xinit'
            xinit = value;
        case 'numiter'
            numIter = value;
        case 'plotrecon'
            plotRecon = value;
        case 'tol'
            tol = value;
        otherwise
            error('pcgWrap: input is not recognized!');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Operator for vectorizing
vecOp = @(x) x(:);

%%% Shape inputs and outputs
AA = @(u) vecOp(A(reshape(u, [Ny, Nx])));
bb = b(:);
xinit = xinit(:);

%%% Run PCG
[utotDom, flag, ~, ~, resvec] = pcg(AA, bb, tol, numIter, [], [], xinit);

utotDom = reshape(utotDom, [Ny, Nx]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outs.flag = flag;
outs.residual = resvec./norm(bb);
iterNum = numel(outs.residual);

if(plotRecon)
    figure(1001);
    set(gcf, 'Color', 'w', 'Name', 'PCG');
    subplot(2, 2, 1);
    imagesc(abs(utotDom));
    axis equal off xy;
    title('abs-utot');
    set(gca, 'FontSize', 16);
    colorbar;
    
    subplot(2, 2, 2);
    imagesc(angle(utotDom));
    axis equal off xy;
    title('angle-utot');
    set(gca, 'FontSize', 16);
    colorbar;
    
    subplot(2, 2, 3:4);
    semilogy(1:iterNum, outs.residual(1:iterNum), 'b-', 'LineWidth', 1.5);
    xlim([1 iterNum]);
    grid on;
    title(sprintf('[Residual: %.2e][Flag: %d]',...
        outs.residual(iterNum), outs.flag));
    set(gca, 'FontSize', 16);
    
    drawnow;
end