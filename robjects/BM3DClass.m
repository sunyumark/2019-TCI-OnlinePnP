classdef BM3DClass < RegularizerClass
    % Class for implementing non-negative constraint.
    %
    % U. S. Kamilov, CIG, WUSTL, 2017.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % public properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        sigSize; % size of the signal
        tau; % regularization parameter
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        function this = BM3DClass(sigSize, tau)
            %  constuctor
            
            this.sigSize = sigSize;
            this.tau = tau;
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
            
            [~, x] = BM3D(1, z, 3e2*this.tau*ratio, 'np', 1);
            
            % do not touch the memory variable
            pout = pin;
        end
        
    end
    
end