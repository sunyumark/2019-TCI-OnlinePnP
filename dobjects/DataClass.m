classdef DataClass < handle
    % DataClass abstract base class for specifying data-fidelity
    % objects for an image reconstruction algorithm.  Data-fidelity
    % classes are specified by deriving from this base class and
    % implementing the methods below.
    %
    % U. S. Kamilov, CIG, WUSTL, 2017.
    
    methods(Abstract)
        % Size of the signal
        sigSize = size(this)
        
        % Evaluate data-fidelity objective
        d = eval(this, x)
        
        % Evaluate gradient
        [g, d] = grad(this, x)
        
        % Draw the signal
        draw(this, x)
    end
end