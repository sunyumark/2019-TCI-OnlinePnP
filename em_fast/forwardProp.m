function utotDomSet = forwardProp(uincDomSet, f, domainGreensFunctionSet,...
    utotDomSet0, dx, dy)
%%% Solves A*utot = uinc for all transmissions and frequencies
%%%
%%% U. S. Kamilov, MERL, 2017.

%%% Dimensions
[~, ~, numTrans, numFreq] = size(uincDomSet);

%%% Initialize the total field
utotDomSet = utotDomSet0;

%%% Convergence flags
convFlags = zeros(numFreq, numTrans);

%%% Loop through frequencies and transmissions
for indFreq = 1:numFreq
    
    domainGreensFunction = domainGreensFunctionSet(:,:,indFreq);
    
    for indTrans = 1:numTrans
        
        %%% Extract incident field
        uincDom = uincDomSet(:,:,indTrans,indFreq);
        utotDom0 = utotDomSet0(:,:,indTrans,indFreq);
        
        %%% Scattering operators
        A = @(u) A_forward(u, f, domainGreensFunction, dx, dy);
        AT = @(z) A_adjoint(z, f, domainGreensFunction, dx, dy);
        
        %%% Compute total field
        [utotDom, outs] = pcgWrap(@(u) AT(A(u)), AT(uincDom),...
            'xinit', utotDom0);
        
        %%% Store convergence flags
        convFlags(indFreq, indTrans) = outs.flag;
        
        %%% Store
        utotDomSet(:,:,indTrans,indFreq) = utotDom;
        
    end
end