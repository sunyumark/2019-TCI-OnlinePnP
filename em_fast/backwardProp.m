function uincDomSet = backwardProp(utotDomSet, f, domainGreensFunctionSet,...
    uincDomSet0, dx, dy)
%%% Solves AH*uinc = utot for all transmissions and frequencies
%%%
%%% U. S. Kamilov, MERL, 2017.

%%% Dimensions
[~, ~, numTrans, numFreq] = size(utotDomSet);

%%% Initialize the total field
uincDomSet = uincDomSet0;

%%% Convergence flags
convFlags = zeros(numFreq, numTrans);

%%% Loop through frequencies and transmissions
for indFreq = 1:numFreq
    
    domainGreensFunction = domainGreensFunctionSet(:,:,indFreq);
    
    for indTrans = 1:numTrans
        
        %%% Extract incident field
        utotDom = utotDomSet(:,:,indTrans,indFreq);
        uincDom0 = uincDomSet0(:,:,indTrans,indFreq);
        
        %%% Scattering operators
        A = @(u) A_forward(u, f, domainGreensFunction, dx, dy);
        AT = @(z) A_adjoint(z, f, domainGreensFunction, dx, dy);
        
        %%% Compute total field
        [uincDom, outs] = pcgWrap(@(u) A(AT(u)), A(utotDom),...
            'xinit', uincDom0);
        
        %%% Store convergence flags
        convFlags(indFreq, indTrans) = outs.flag;
        
        %%% Store
        uincDomSet(:,:,indTrans,indFreq) = uincDom;
    end
end
