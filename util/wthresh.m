function xhat = wthresh(y, tau)
%%% Soft thresholding function
%%% U. S. Kamilov, 2015.

% compute the norm in dimension 3
norm_y = abs(y);

amp = max(norm_y-tau, 0);
norm_y(norm_y<=0) = 1; % to avoid division by 0

xhat = amp./norm_y .* y;