function z = H_forward(x, g, dx, dy)
%%% Computes the scattered field at the sensor locations specified by the
%%% set of Green's functions.
%%%
%%% Input:
%%% - x: (Ny x Nx) contrast-source
%%% - g: (Ny x Nx x Nr) Green's functions
%%% - dx, dy: sampling steps
%%%
%%% Ouput:
%%% - z: (Nr x 1) scattered field
%%%
%%% U. S. Kamilov, MERL, 2017.

%%% number of transmissions
Nr = size(g, 3);

%%% contrast-source
x = repmat(x,[1,1,Nr]);

%%% compute propagated field
z = dx*dy*sum(sum(g.*x,1),2); % measure
z = z(:);