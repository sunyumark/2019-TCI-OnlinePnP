%%% Script for testing the 2D convolution
%%%
%%% U. S. Kamilov, MERL, 2017.

close all; clc;

%%

Nx = 16;
Ny = 16;

f = randn(Ny, Nx);
g = randn(2*size(f));
Hf = conv2D(f, g);

z = randn(size(Hf));
HTz = conv2DAdj(z, g);

z(:)' * Hf(:)
HTz(:)' * f(:)

%%

Nx = 16;
Ny = 16;
dx = 0.1;
dy = 0.1;

f = randn(Ny, Nx)+1j*randn(Ny, Nx);
g = randn(2*size(f))+1j*randn(2*size(f));
u = randn(Ny, Nx)+1j*randn(Ny, Nx);

Au = A_forward(u, f, g, dx, dy);

z = randn(size(Au))+1j*randn(size(Au));

AHz = A_adjoint(z, f, g, dx, dy);

z(:)' * Au(:)
AHz(:)' * u(:)

%%

Nx = 16;
Ny = 16;
Nt = 8;
Nr = 8;
Nf = 16;
dx = 0.1;
dy = 0.1;

x = randn(Ny, Nx, Nt, Nf)+1j*randn(Ny, Nx, Nt, Nf);
g = randn(2*Ny, 2*Nx, 2*Nf)+1j*randn(2*Ny, 2*Nx, 2*Nf);

Gx = G_forward_full(x, g, dx, dy);

z = randn(size(Gx))+1j*randn(size(Gx));

GTz = G_adjoint_full(z, g, dx, dy);

z(:)' * Gx(:)
GTz(:)' * x(:)

%%
Nx = 16;
Ny = 16;
Nr = 8;
dx = 0.1;
dy = 0.1;

f = randn(Ny, Nx)+1j*randn(Ny, Nx);
g = randn(Ny, Nx, Nr)+1j*randn(Ny, Nx, Nr);
uin = randn(Ny, Nx)+1j*randn(Ny, Nx);

Hf = propagateToSensor(f, uin, g, dx, dy);

z = randn(size(Hf))+1j*randn(size(Hf));
HTz = propagateToSensorAdj(z, uin, g, dx, dy);

z(:)' * Hf(:)
HTz(:)' * f(:)

%%
Nx = 16;
Ny = 16;
Nr = 8;
dx = 0.1;
dy = 0.1;

x = randn(Ny, Nx)+1j*randn(Ny, Nx);
g = randn(Ny, Nx, Nr)+1j*randn(Ny, Nx, Nr);

Hx = H_forward(x, g, dx, dy);

z = randn(size(Hx))+1j*randn(size(Hx));
HTz = H_adjoint(z, g, dx, dy);

z(:)' * Hx(:)
HTz(:)' * x(:)

%%
Nx = 16;
Ny = 16;
Nt = 8;
Nr = 8;
Nf = 16;
dx = 0.1;
dy = 0.1;

x = randn(Ny, Nx, Nt, Nf)+1j*randn(Ny, Nx, Nt, Nf);
g = randn(Ny, Nx, Nr, Nf)+1j*randn(Ny, Nx, Nr, Nf);
receiverMaskSet = randi(2, [Nt, Nr, Nf])-1;

Hx = H_forward_full(x, g, receiverMaskSet, dx, dy);

z = randn(size(Hx))+1j*randn(size(Hx));
HTz = H_adjoint_full(z, g, receiverMaskSet, dx, dy);

z(:)' * Hx(:)
HTz(:)' * x(:)

%% FullPropagate to Sensor

Nx = 16;
Ny = 16;
Nr = 8;
Nt = 4;
Nf = 4;

f = randn(Ny, Nx)+1j*randn(Ny, Nx);
uin = randn(Ny, Nx, Nt, Nf)+1j*randn(Ny, Nx, Nt, Nf);
g = randn(Ny, Nx, Nr, Nf)+1j*randn(Ny, Nx, Nr, Nf);
receiverMaskSet = randi(2, [Nt, Nr, Nf])-1;

Hf = fullPropagateToSensor(f, uin, g, receiverMaskSet, dx, dy);

z = randn(size(Hf))+1j*randn(size(Hf));
HTz = fullPropagateToSensorAdj(z, uin, g, receiverMaskSet, dx, dy);

z(:)' * Hf(:)
HTz(:)' * f(:)

%% FullPropagate to Sensor

Nx = 16;
Ny = 16;
dx = 0.1;
dy = 0.1;

f = randn(Ny, Nx)+1j*randn(Ny, Nx);
g = randn(2*size(f))+1j*randn(2*size(f));
x = randn(Ny, Nx)+1j*randn(Ny, Nx);

Bx = forwardProp(x, f, g, dx, dy);

z = randn(size(Bx))+1j*randn(size(Bx));
BTz = backwardProp(z, f, g, dx, dy);

z(:)' * Bx(:)
BTz(:)' * x(:)


