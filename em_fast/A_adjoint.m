function u = A_adjoint(z, f, g, dx, dy)
%%% Adjoint of A = I - G*diag(f) operator required for forward scattering
%%%
%%% Input:
%%% - z: (Ny x Nx) field
%%% - f: (Ny x Nx) object
%%% - g: 2*(Ny x Nx) Green's function
%%%
%%% Output:
%%% - u: (Ny x Nx) field
%%%
%%% U. S. Kamilov, MERL, 2017.

assert(all(size(g)==2*size(z)), 'all(size(g)==2*size(z))');
assert(all(size(f)==size(z)), 'all(size(f)==size(z))');

u = z - dx*dy*conj(f).*conv2DAdj(z, g);