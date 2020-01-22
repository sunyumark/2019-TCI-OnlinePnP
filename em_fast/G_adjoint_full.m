function xSet = G_adjoint_full(zSet, domainGreensFunctionSet, dx, dy)
%%% Convolve with adjoint domain Green's function

%%% Dimensions
numTrans = size(zSet, 3);
numFreq = size(zSet, 4);

%%% Initizalize
xSet = zeros(size(zSet));

for indFreq = 1:numFreq
    
    g = domainGreensFunctionSet(:,:,indFreq);
    
    for indTrans = 1:numTrans
        
        %%% Extract
        z = zSet(:,:,indTrans,indFreq);
        
        %%% Predicted scattered field
        x = dx*dy*conv2DAdj(z, g);
        
        %%% Store
        xSet(:,:,indTrans,indFreq) = x;
        
    end
end