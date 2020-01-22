function zSet = H_forward_full(xSet, sensorGreensFunctionSet, receiverMaskSet, dx, dy)
%%% Computes the scattered field at the sensor locations specified by the
%%% set of Green's functions.
%%% U. S. Kamilov, MERL, 2017.

%%% Dimensions
numTrans = size(xSet, 3);
numFreq = size(xSet, 4);
numRec = size(sensorGreensFunctionSet, 3);

%%% Initialize
zSet = zeros(numTrans, numRec, numFreq);

for indFreq = 1:numFreq
    
    g = sensorGreensFunctionSet(:,:,:,indFreq);
    
    for indTrans = 1:numTrans
        
        %%% Extract
        receiverMask = receiverMaskSet(indTrans,:,indFreq);
        x = xSet(:,:,indTrans,indFreq);
        
        %%% Predicted scattered field
        z = H_forward(x, g, dx, dy);
        z = receiverMask(:) .* z;
        
        %%% Store
        zSet(indTrans,:,indFreq) = z;
        
    end
end