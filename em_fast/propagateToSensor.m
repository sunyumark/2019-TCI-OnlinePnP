function uscat = propagateToSensor(f, uin, g, dx, dy)
%%% Computes the scattered field at the sensor locations specified by the
%%% set of Green's functions.
%%%
%%% Input:
%%% - f: (Ny x Nx) scattering potential
%%% - uin: (Ny x Nx x Nt) input field
%%% - g: (Ny x Nx x Nr) Green's functions
%%% - dx, dy: sampling steps
%%%
%%% Ouput:
%%% - z: (Nr x 1) scattered field
%%%
%%% U. S. Kamilov, MERL, 2017.

%%% number of transmissions
Nr = size(g, 3);
[Ny,Nx,Nt] = size(uin);

%%% reshape matrix
f = repmat(f,[1,1,Nt]);
f = reshape(f,[Ny*Nx,Nt]);
uin = reshape(uin,[Ny*Nx,Nt]);
g = reshape(g,[Ny*Nx,Nr]);

%%% contrast-source
contSrc = f .* uin;  % [Ny*Nx, Nt]

%%% compute propagated field
uscat = dx*dy*conj(contSrc')*g; % measure [Nt, Nr]
