function f = conv2DAdj(z, g)
%%% Convolves with the adjoint of the function g.
%%%
%%% Input:
%%% - z: (Ny x Nx) object
%%% - g: 2*(Ny x Nx) function inside the object
%%%
%%% Output:
%%% - f: (Ny x Nx) output
%%%
%%% U. S. Kamilov, MERL, 2017.

%%% Validate the size
assert(all((2*size(z))==size(g)), 'all((2*size(z))==size(g))');

%%% Shift
g = circshift(g, -size(g)/2);

%%% size of computational domain and number of inputs
[Ny, Nx] = size(z);

%%% zero post-padding
z = padarray(z, [Ny, Nx], 0, 'post');

%%% convolve
f = ifft2(fft2(z).*conj(fft2(g)));

%%% remove zero-padding
f = f(1:Ny,1:Nx);