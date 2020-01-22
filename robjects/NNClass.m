classdef NNClass < RegularizerClass
    % Class for implementing non-negative constraint.
    %
    % U. S. Kamilov, CIG, WUSTL, 2017.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % public properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        sigSize; % size of the signal
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        function this = NNClass(sigSize)
            %  constuctor
            
            this.sigSize = sigSize;
        end
        
        function p = init(this)
            % does nothing
            
            p = zeros(this.sigSize);
        end
        
        function r = eval(~, x)
            % evaluate l1-norm
            
            r = any(x(:) < 0)*1e3;
        end
        
        function [x, pout] = prox(~, z, ~, pin)
            % soft-thresholding
            
            % project
            x = z;
            x(x<0) = 0;
            
            % do not touch the memory variable
            pout = pin;
        end
        
    end
    
end