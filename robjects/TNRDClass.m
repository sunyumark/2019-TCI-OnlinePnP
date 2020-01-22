classdef TNRDClass < RegularizerClass
    % Class for implementing non-negative constraint.
    %
    % Yu Sun, CIG, WUSTL, 2017.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % public properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        sigSize; % size of the signal
        tnrd;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        function this = TNRDClass(sigSize, path)
            %  constuctor
            
            this.sigSize = sigSize;
            this.tnrd = load(path);
            
        end
        
        function p = init(this)
            % does nothing
            
            p = zeros(this.sigSize);
        end
        
        function r = eval(~, x)
            % evaluate l1-norm
            
            r = 0;
        end
        
        function [x, pout] = prox(this, z, step, pin, ratio)
            % soft-thresholding
            
            % project
            z(z<0) = 0;
            
            x = TNRD(z, this.tnrd);
            
            % do not touch the memory variable
            pout = pin;
        end
        
    end
    
end