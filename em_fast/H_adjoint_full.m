function xSet = H_adjoint_full(zSet, sensorGreensFunctionSet, receiverMaskSet, dx, dy)
%%% Propagate all the scattered fields to computational domain
%%%
%%% U. S. Kamilov, MERL, 2017.

%%% Dimensions
numTrans = size(zSet, 1);
numFreq = size(zSet, 3);
Ny = size(sensorGreensFunctionSet, 1);
Nx = size(sensorGreensFunctionSet, 2);

%%% Initialize
xSet = zeros(Ny, Nx, numTrans, numFreq);

for indFreq = 1:numFreq
    
    g = sensorGreensFunctionSet(:,:,:,indFreq);
    
    for indTrans = 1:numTrans
        
        %%% Extract
        receiverMask = receiverMaskSet(indTrans,:,indFreq);
        z = zSet(indTrans, :, indFreq);
        z = receiverMask(:) .* z(:);
        
        %%% Predicted scattered field
        x = H_adjoint(z, g, dx, dy);
        
        %%% Add contribution
        xSet(:,:,indTrans,indFreq) = x;
    end
end