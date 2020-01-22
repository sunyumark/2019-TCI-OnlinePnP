function z = conv2D(f, g)
%%% Convolves with the function g.
%%%
%%% Input:
%%% - f: (Ny x Nx) object
%%% - g: 2*(Ny x Nx) function inside the object
%%%
%%% Output:
%%% - z: (Ny x Nx) output
%%%
%%% U. S. Kamilov, MERL, 2017.

%%% Validate the size
assert(all((2*size(f))==size(g)), 'all((2*size(f))==size(g))');

%%% Shift
g = circshift(g, -size(g)/2);

%%% size of computational domain and number of inputs
[Ny, Nx] = size(f);

%%% zero post-padding
f = padarray(f, [Ny, Nx], 0, 'post');

%%% convolve
z = ifft2(fft2(f).*fft2(g));

%%% remove zero-padding
z = z(1:Ny,1:Nx);