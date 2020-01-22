classdef FirstBornClass < DataClass
    % Class for wave scattering based on the first Born approximation.
    %
    % U. S. Kamilov, CIG, WUSTL, 2018.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % public properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        y; % measurements
        sigSize; % size of the signal
        emParams; % scattering parameters
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        function this = FirstBornClass(y, emParams)
            % constructor
            this.y = y;
            this.emParams = emParams;
            this.sigSize = [emParams.Ny, emParams.Nx];
        end
        
        function sigSize = size(this)
            % size of the image
            
            sigSize = this.sigSize;
        end
        
        function d = eval(this, x)
            % evaluate data-fidelity cost
            
            z = this.fmult(x, this.emParams);
            d = 0.5*norm(this.y(:)-z(:))^2;
        end
        
        function [g] = grad(this, x)
            % evaluate gradient and data-fidelity
            
            z = this.fmult(x, this.emParams);
            g = this.ftran(z-this.y, this.emParams);
            g = real(g); % only keep the real part
%             d = 0.5*norm(this.y(:)-z(:))^2;
        end
        
        function [gSub, d, keepIdx] = gradSubFull(this, x, keepIdxPrev)
            % evaluate stochastic gradient and data-fidelity
            % get randomly sub-sampled em and keepIdx
            [emSub, keepIdx] = this.stochastizeEm(this.emParams, keepIdxPrev);
            % drop y
            ySub = this.y(keepIdx,:);
            % drop z
            z = this.fmult(x, this.emParams);
            zSub = z(keepIdx,:);
            % drop 
            gSub = this.ftran(zSub-ySub, emSub);
            gSub = real(gSub); % only keep the real part
            d = 0.5*norm(this.y(:)-z(:))^2;
        end
        
        function [gSub, dSub, keepIdx] = gradSub(this, x, keepIdxPrev)
            % evaluate stochastic gradient
            % get randomly sub-sampled em and keepIdx
            [emSub, keepIdx] = this.stochastizeEm(this.emParams, keepIdxPrev);
            % drop y
            ySub = this.y(keepIdx,:);
            % drop z
            zSub = this.fmult(x, emSub);
            % drop 
            gSub = this.ftran(zSub-ySub, emSub);
            gSub = real(gSub); % only keep the real part
            dSub = 0.5*norm(ySub(:)-zSub(:))^2;
        end
        
        function draw(~, x)
            % show real(x) as an image            
            imagesc(real(x));
            axis equal off tight;
            colormap gray;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % static methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods(Static)
        
        function [emParamsOut, keepIdx] = stochastizeEm(emParams, keepIdxPrev)
            % randomly sub-samples the measurements
            % related variables
            keepIdx = mod(keepIdxPrev, 60)+1;
            emParamsOut = emParams;
            emParamsOut.numTrans = length(keepIdx);
            emParamsOut.uincDom = emParams.uincDom(:,:,keepIdx);
            emParamsOut.receiverMask = emParams.receiverMask(keepIdx,:);
        end
        
        function z = fmult(x, emParams)
            % forward operator
            %
            % Input:
            % - x: 2D image
            % - emParams: scattering parameters
            %
            % Output:
            % - z: scattered field at each sensor
            %
            % U. S. Kamilov, CIG, WUSTL, 2018.
            
            z = fullPropagateToSensor(x,...
                emParams.uincDom, emParams.sensorGreensFunction,...
                emParams.receiverMask, emParams.dx, emParams.dy);
            
        end
        
        function x = ftran(z, emParams)
            % adjoint operator
            %
            % Input:
            % - z: scattered field at each sensor
            % - emParams: scattering parameters
            %
            % Output:
            % - x: 2D image
            %
            % U. S. Kamilov, CIG, WUSTL, 2018.
            
            x = fullPropagateToSensorAdj(z,...
                emParams.uincDom, emParams.sensorGreensFunction,...
                emParams.receiverMask, emParams.dx, emParams.dy);
            
        end
    end
end