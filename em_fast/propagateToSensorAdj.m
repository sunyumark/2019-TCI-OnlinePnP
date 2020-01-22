function fhat = propagateToSensorAdj(uscat, uin, g, dx, dy)
%%% Adjoint of the propagation function
%%%
%%% Inputs:
%%% - uscat: (Nt x Nr) scattered field
%%% - uin: (Ny x Nx x Nt) input field inside computational domain
%%% - g: (Ny x Nx x Nr) Green's function
%%% - dx,dy: sampling steps
%%%
%%% Outputs:
%%% - fhat: (Ny x Nx x Nt) scattering potential
%%%
%%% U. S. Kamilov, MERL, 2016

%%% dimensions of the problem
Nr = size(g,3);
[Ny,Nx,Nt] = size(uin);

%%% reorganize scattered field
g = reshape(g,[Ny*Nx,Nr]);
uin = reshape(uin,[Ny*Nx,Nt]);

%%% compute inner product
contSrc = dx*dy*conj(g)*conj(uscat');
contSrc = reshape(contSrc, [Ny,Nx,Nt]);
uin = reshape(uin, [Ny,Nx,Nt]);

%%% Multiply to adjoint field
fhat = conj(uin).*contSrc;