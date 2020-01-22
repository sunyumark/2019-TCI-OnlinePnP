function x = H_adjoint(z, g, dx, dy)
%%% Adjoint of the propagation function
%%%
%%% Inputs:
%%% - z: (Nr x 1) scattered field
%%% - g: (Ny x Nx x Nr) Green's function
%%% - dx,dy: sampling steps
%%%
%%% Outputs:
%%% - x: (Ny x Nx) scattering potential
%%%
%%% U. S. Kamilov, MERL, 2016

%%% dimensions of the problem
Nr = size(g,3);
Ny = size(g,1);
Nx = size(g,2);

%%% reorganize scattered field
z = reshape(z,[1,1,Nr]);
z = repmat(z,[Ny,Nx]);

%%% compute inner product
x = dx*dy*sum(z.*conj(g),3);