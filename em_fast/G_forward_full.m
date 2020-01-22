function zSet = G_forward_full(xSet, domainGreensFunctionSet, dx, dy)
%%% Convolve with domain Green's function

%%% Dimensions
numTrans = size(xSet, 3);
numFreq = size(xSet, 4);

%%% Initizalize
zSet = zeros(size(xSet));

for indFreq = 1:numFreq
    
    g = domainGreensFunctionSet(:,:,indFreq);
    
    for indTrans = 1:numTrans
        
        %%% Extract
        x = xSet(:,:,indTrans,indFreq);
        
        %%% Predicted scattered field
        z = dx*dy*conv2D(x, g);
        
        %%% Store
        zSet(:,:,indTrans,indFreq) = z;
        
    end
end