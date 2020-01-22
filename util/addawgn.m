function [y, noise] = addawgn(x, inputSnr)
% ADDAWGN returns an array
%   y = x + noise
% where noise is an additive white Gaussian noise of variance
% noiseVariance.
% The formula for computation is
%   10*log10([norm(x(:))^2]/[N*noiseVariance]) = ISNR
%
% U. S. Kamilov, MERL, 2016

noiseNorm = norm(x(:)) * 10^(-inputSnr/20); % compute the desired norm

% generate AWGN
if(isreal(x))
    noise = randn(size(x)); 
else
    noise = randn(size(x)) + 1j*randn(size(x));
end


noise = noise/norm(noise(:)) * noiseNorm; % adjust the noise power


y = x + noise; % add noise