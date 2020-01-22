function tau = optimizeTau(x, algoHandle, taurange)
%%% Optimize regularization parameter

evaluateSNR = @(x, xhat) 20*log10(norm(x(:))/norm(x(:)-xhat(:)));

options = optimset('Display', 'iter', 'TolX', 1e-16, 'MaxFunEvals', 20);

fun = @(tau) -evaluateSNR(x, algoHandle(tau));

tau = fminbnd(fun, taurange(1), taurange(2), options);