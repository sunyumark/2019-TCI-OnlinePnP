classdef CisorClass < DataClass
    % Class for wave scattering based on CISOR.
    %
    % U. S. Kamilov, CIG, WUSTL, 2018.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % public properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        y; % measurements
        sigSize; % size of the signal
        emParams; % scattering parameters
        utot; % total field u
        vtot; % total field v
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % public methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        function this = CisorClass(y, emParams)
            % constructor
            this.y = y;
            this.emParams = emParams;
            this.sigSize = [emParams.Ny, emParams.Nx];
            this.utot = zeros(size(emParams.uincDom));
            this.vtot = zeros(size(emParams.uincDom));
        end
        
        function sigSize = size(this)
            % size of the image
            
            sigSize = this.sigSize;
        end
        
        function d = eval(this, x)
            % evaluate data-fidelity cost
            
            [z, this.utot] = this.forward(x, this.utot, this.emParams);
            d = 0.5*norm(this.y(:)-z(:))^2;
        end
        
        function [g, d] = grad(this, x)
            % evaluate the gradient and data-fidelity
            
            % objective
            [z, this.utot] = this.forward(x, this.utot, this.emParams);
            d = 0.5*norm(this.y(:)-z(:))^2;
            
            % gradient
            g = this.backward(z-this.y, x,...
                this.utot, this.vtot, this.emParams);
        end
        
        function [gSub, d] = gradSub(this, x, keepNum)
            % evaluate the gradient and data-fidelity
            
            % objective
            [z, this.utot] = this.forward(x, this.utot, this.emParams);
            d = 0.5*norm(this.y(:)-z(:))^2;
            
            % sub-sample measurements
            [emSub, keepIdx] = this.stochastizeEm(this.emParams, keepNum);
            zSub = z(keepIdx,:);
            ySub = this.y(keepIdx,:);
            utotSub = this.utot(:,:,keepIdx);
            vtotSub = this.vtot(:,:,keepIdx);
            
            % gradient
            gSub = this.backward(zSub-ySub, x,...
                utotSub, vtotSub, emSub);
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
        
        function [emParamsOut, keepIdx] = stochastizeEm(emParams, keepNum)
            % randomly sub-samples the measurements
            % related variables
            keepIdx = randperm(emParams.numTrans, keepNum);
            emParamsOut = emParams;
            emParamsOut.numTrans = keepNum;
            emParamsOut.uincDom = emParams.uincDom(:,:,keepIdx);
            emParamsOut.receiverMask = emParams.receiverMask(keepIdx,:);
        end
        
        function [z, u] = forward(x, u, emParams)
            % given the object compute the scattered field
            
            % extract the incident field
            u0 = emParams.uincDom;
            
            % compute total field with CISOR
            u = forwardProp(u0, x, emParams.domainGreensFunction, u,...
                emParams.dx, emParams.dy);
            
            % propagate to sensor
            z = fullPropagateToSensor(x, u,...
                emParams.sensorGreensFunction,...
                emParams.receiverMask, emParams.dx, emParams.dy);
        end
        
        function g = backward(r, x, u, v, emParams)
            % given the residual compute the gradient
            %  to understand see eq. (12) and (13) in the paper.
            
            HTr = H_adjoint_full(r,...
                emParams.sensorGreensFunction,...
                emParams.receiverMask, emParams.dx, emParams.dy);
            
            xHTr = HTr .* repmat(x, [1, 1, size(HTr, 3)]);
            
            v = backwardProp(xHTr, x, emParams.domainGreensFunction, v,...
                emParams.dx, emParams.dy);
            
            GTv = G_adjoint_full(v, emParams.domainGreensFunction,...
                emParams.dx, emParams.dy);
            
            g = sum(real(conj(u) .* (HTr + GTv)), 3);
        end
        
    end
end