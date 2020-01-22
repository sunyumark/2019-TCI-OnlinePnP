function z = A_forward(u, f, g, dx, dy)
%%% A*u = (I - G*diag(f))*u operator required for forward scattering
%%%
%%% Input:
%%% - u: (Ny x Nx) field
%%% - f: (Ny x Nx) object
%%% - g: 2*(Ny x Nx) Green's function
%%%
%%% Output:
%%% - z: (Ny x Nx) field
%%%
%%% U. S. Kamilov, MERL, 2017.

assert(all(size(g)==2*size(u)), 'all(size(g)==2*size(u))');
assert(all(size(f)==size(u)), 'all(size(f)==size(u))');

z = u - dx*dy*conv2D(f.*u, g);