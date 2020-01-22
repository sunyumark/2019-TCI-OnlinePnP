classdef RegularizerClass < handle
    % RegularizerClass abstract base class for specifying regularizer
    % objects for an image reconstruction algorithm.  Regularizer
    % classes are specified by deriving from this base class and
    % implementing the methods below.
    %
    % U. S. Kamilov, CIG, WUSTL, 2017.
    
    methods(Abstract)
        % Initialize memory variable
        p = init(this)
        
        % Evaluate regularizer objective
        r = eval(this, x)
        
        % Evaluate proximal
        [x, pout] = prox(this, z, step, pin)
    end
    
end