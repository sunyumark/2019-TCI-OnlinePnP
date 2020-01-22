function lam = powerIter(A, sigSize, varargin)
% POWERITER finds the largest eigen-value of the matrix
%   A = H^T * H
%
% Input:
%   - A: forward model operator
%   - sigSize: size of the vectors
%   - varargin: options for the algorithm
%       Options:
%       - numIter: number of iterations (def: 100)
%       - tol: stopping criterion (def: 1e-6)
%       - verbose: print command line message (def: false)
%
% Output:
%   - lam: lamda_max(HTH)
%
% U. S. Kamilov, MERL, 2015.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numIter = 100;
tol = 1e-6;
verbose = false;

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
        case 'tol'
            tol = value;
        case 'verbose'
            verbose = value;
        otherwise
            error('powerIter: input is not recognized!');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shortcuts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Shortcuts
x = randn(sigSize);
x = x/norm(x(:));

lam = 1;

for indIter = 1:numIter
    xnext = A(x);
    
    lamNext = real((xnext(:)' * x(:))/(x(:)' * x(:)));
    
    xnext = xnext/norm(xnext(:));
    
    relDiff = abs(lamNext-lam)/abs(lam);
    
    x = xnext;
    lam = lamNext;
    
    if(verbose)
        fprintf('[%d/%d] lam = %.6e, relDiff = %.1e\n',...
            indIter, numIter, lam, relDiff);
    end
    
    if(relDiff < tol)
        break;
    end
end