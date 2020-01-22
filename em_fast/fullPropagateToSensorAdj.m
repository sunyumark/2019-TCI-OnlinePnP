function f = fullPropagateToSensorAdj(uscatPredSet, utotDomSet,...
    sensorGreensFunctionSet, receiverMaskSet, dx, dy)
%%% Propagate all the scattered fields to computational domain
%%%
%%% U. S. Kamilov, MERL, 2017.

%%% Dimensions
numFreq = size(uscatPredSet, 3);
Ny = size(sensorGreensFunctionSet, 1);
Nx = size(sensorGreensFunctionSet, 2);

%%% Initialize
f = zeros(Ny, Nx);

for indFreq = 1:numFreq
    
    sensorGreensFunction = sensorGreensFunctionSet(:,:,:,indFreq);
        
    %%% Extract
    receiverMask = receiverMaskSet(:,:,indFreq);
    utotDom = utotDomSet(:,:,:,indFreq);
    uscatPred = uscatPredSet(:,:,indFreq);
        
    %%% Predicted scattered field
    uscatPred = uscatPred .* receiverMask;
    contSrc = propagateToSensorAdj(uscatPred, utotDom, sensorGreensFunction, dx, dy);
        
    %%% Save each contribution
    f(:,:) = sum(contSrc,3);

end