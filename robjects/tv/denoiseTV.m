function [x,P,iter,L]=denoiseTV(y,lambda,varargin)

%This method uses the TV regularizer and can be applied only to 2D data.
%x: denoised image
%P: dual variables
%L: Lipschitz constant
%iter: number of iterations for getting to the solution.
%bc: Boundary conditions for the differential operators
%    ('reflexive'|'circular'|'zero')

[maxiter,L,tol,optim,verbose,img,bounds,P,bc]=process_options(varargin,...
  'maxiter',100,'L',8,'tol',1e-5,'optim','fgp','verbose',false,'img',[],...
  'bounds',[-inf +inf],'P',zeros([size(y) 2]),'bc','reflexive');

if isempty(L)
  L=Lipschitz(y)/1.25;% Lipschitz constant for the TV Operator
end

count=0;
flag=false;

if verbose
  fprintf('******************************************\n');
  fprintf('**     Denoising with TV Regularizer    **\n');
  fprintf('******************************************\n');
  fprintf('#iter     relative-dif   \t fun_val         Duality Gap    \t   ISNR\n')
  fprintf('====================================================================\n');
end
switch optim
  case 'fgp'
    t=1;
    F=P;
    for i=1:maxiter
      K=y-lambda*AdjTVOp2D(F,bc);      
      Pnew=F+(1/(L*lambda))*TVOp2D(project(K,bounds),bc);
      Pnew=projectL2(Pnew);
      
      %relative error
      re=norm(Pnew(:)-P(:))/norm(Pnew(:));
      if (re<tol)
        count=count+1;
      else
        count=0;
      end
      
      tnew=(1+sqrt(1+4*t^2))/2;
      F=Pnew+(t-1)/tnew*(Pnew-P);
      P=Pnew;
      t=tnew;
      
      if verbose
        if ~isempty(img)
          k=y-lambda*AdjTVOp2D(P,bc);
          x=project(k,bounds);
          fun_val=cost(y,x,lambda,bc);
          dual_fun_val=dualcost(y,k,bounds);
          dual_gap=(fun_val-dual_fun_val);
          ISNR=20*log10(norm(y-img,'fro')/norm(x-img,'fro'));
          % printing the information of the current iteration
          fprintf('%3d \t %10.5f \t %10.5f \t %2.8f \t %2.8f\n',i,re,fun_val,dual_gap,ISNR);
        else
          k=y-lambda*AdjTVOp2D(P,bc);
          x=project(k,bounds);
          fun_val=cost(y,x,lambda,bc);
          dual_fun_val=dualcost(y,k,bounds);
          dual_gap=(fun_val-dual_fun_val);
          fprintf('%3d \t %10.5f \t %10.5f \t %2.8f\n',i,re,fun_val,dual_gap);
        end
      end
      
      if count >=5
        flag=true;
        iter=i;
        break;
      end
    end
    
  case 'gp'
    
    for i=1:maxiter
      
      K=y-lambda*AdjTVOp2D(P,bc);
      
      Pnew=P+(1/(L*lambda))*TVOp2D(project(K,bounds),bc);
      Pnew=projectL2(Pnew);
      
      %relative error
      re=norm(Pnew(:)-P(:))/norm(Pnew(:));
      if (re<tol)
        count=count+1;
      else
        count=0;
      end
      
      P=Pnew;
      
      if verbose
        if ~isempty(img)
          k=y-lambda*AdjTVOp2D(P,bc);
          x=project(k,bounds);
          fun_val=cost(y,x,lambda,bc);
          dual_fun_val=dualcost(y,k,bounds);
          dual_gap=(fun_val-dual_fun_val);
          ISNR=20*log10(norm(y-img,'fro')/norm(x-img,'fro'));
          % printing the information of the current iteration
          fprintf('%3d \t %10.5f \t %10.5f \t %2.8f \t %2.8f\n',i,re,fun_val,dual_gap,ISNR);
        else
          k=y-lambda*AdjTVOp2D(P,bc);
          x=project(k,bounds);
          fun_val=cost(y,x,lambda,bc);
          dual_fun_val=dualcost(y,k,bounds);
          dual_gap=(fun_val-dual_fun_val);
          fprintf('%3d \t %10.5f \t %10.5f \t %2.8f\n',i,re,fun_val,dual_gap);
        end
      end
      
      if count >=5
        flag=true;
        iter=i;
        break;
      end
    end
end

if ~flag
  iter=maxiter;
end

x=project(y-lambda*AdjTVOp2D(P,bc),bounds);


function Df=TVOp2D(f,bc) %TV operator with reflexive boundary conditions

[r,c]=size(f);
Df=zeros(r,c,2);
Df(:,:,1)=shift(f,[-1,0],bc)-f;
Df(:,:,2)=shift(f,[0,-1],bc)-f;


function g=AdjTVOp2D(P,bc) %Adjoint TV operator

P1=P(:,:,1);
P1=shiftAdj(P1,[-1,0],bc)-P1;
P2=P(:,:,2);
P2=shiftAdj(P2,[0,-1],bc)-P2;
g=P1+P2;

function PB=projectL2(B)

%Check which cubes (matrices NxMx2) don't belong to the l2-norm ball.
%For those cubes normalize their values by
%PB=B(:,:,1)/max(1,sqrt(B(:,:,1)^2+B(:,:,2)^2),
%B(:,:,2)/max(1,sqrt(B(:,:,1)^2+B(:,:,2)^2)

%K=max(1,sqrt(sum(B.^2,3)));
%PB=B./repmat(K,[1 1 2]);
PB=B./repmat(max(1,sqrt(sum(B.^2,3))),[1 1 2]);

function Px=project(x,bounds)
lb=bounds(1);%lower box bound
ub=bounds(2);%upper box bound

if isequal(lb,-Inf) && isequal(ub,Inf)
  Px=x;
elseif isequal(lb,-Inf) && isfinite(ub)
  x(x>ub)=ub;
  Px=x;
elseif isequal(ub,Inf) && isfinite(lb)
  x(x<lb)=lb;
  Px=x;
else
  x(x<lb)=lb;
  x(x>ub)=ub;
  Px=x;
end

function [Q,TVnorm]=cost(y,f,lambda,bc)

fx=shift(f,[-1,0],bc)-f;
fy=shift(f,[0,-1],bc)-f;

TVf=sqrt(fx.^2+fy.^2);% Amplitude of the gradient vector

TVnorm=sum(TVf(:));
Q=0.5*norm(y-f,'fro')^2+lambda*TVnorm;

function Q=dualcost(y,f,bounds)
r=f-project(f,bounds);
Q=0.5*(sum(r(:).^2)+sum(y(:).^2)-sum(f(:).^2));


function L=Lipschitz(y)

%Finds the Lipschitz constant for the function f(A)=0.5*||L*(A)-y||^2,
%where L* is the adjoint Laplacian operator and A is an image, e.g. H(A),
%where H is the Laplacian operator.

[r,c]=size(y);
%The Lipschitz constant for f(A) is equal to ||LL*||=||L||^2
%where r is the spectral radius of the operator.
%LL*=(Dxx+Dyy)*(Dxx+Dyy)*, which is a circulant operator since each one
%of the suboperators are circulant. The addition and multiplication of
%circulant operators results to a circulant operator.

hx=[1 -1 0 ]';% Dx Operator
hx=padarray(hx,[0 1]);
hy=hx'; %Dy Operator
center=[2 2];

hx(r,c)=0;hy(r,c)=0;

%Circular shift of the operator T1 so as to have consistent results by
%applying either one of the following 2 operations on an image x (T1zp is
%the zero padded operator)
%1)imfilter(x,T1,'conv','circular')
%2)ifft2(fft2(x).*fft2(T1zp));
hx=circshift(hx,1-center);
hy=circshift(hy,1-center);


%Operator eigenvalues
Op_eig=abs(fft2(hx)).^2+abs(fft2(hy)).^2;
L=max(Op_eig(:));
