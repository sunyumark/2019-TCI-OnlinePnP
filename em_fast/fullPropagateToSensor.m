function uscatPredSet = fullPropagateToSensor(f, utotDomSet,...
    sensorGreensFunctionSet, receiverMaskSet, dx, dy)
%%% Propagate all the total fields to the sensors
%%%
%%% U. S. Kamilov, MERL, 2017.

%%% Dimensions
numTrans = size(utotDomSet, 3);
numFreq = size(utotDomSet, 4);
numRec = size(sensorGreensFunctionSet, 3);

%%% Initialize
uscatPredSet = zeros(numTrans, numRec, numFreq);

for indFreq = 1:numFreq
    
    sensorGreensFunction = sensorGreensFunctionSet(:,:,:,indFreq);
        
    receiverMask = receiverMaskSet(:,:,indFreq); % [Nr, Nt, Nf]
    utotDom = utotDomSet(:,:,:,indFreq);
        
    %%% Predicted scattered field
    uscatPred = propagateToSensor(f, utotDom, sensorGreensFunction, dx, dy);
    uscatPred = receiverMask .* uscatPred;
        
    %%% Store
    uscatPredSet(:,:,indFreq) = uscatPred;

end