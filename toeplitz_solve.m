function x = toeplitz_solve(c, r, b)

n = length(c);

omega = exp(1i * pi / n .* (0 : n - 1));

G = [ 1 , 2*r(1) ; zeros(n-1, 1), r(end:-1:2).' + c(2:end) ]; % G(1,2) = 2 * c(1);
H = [ c(end:-1:2) - r(2:end).' , zeros(n-1,1) ; 0 1 ];

Gh = ifft(G) * sqrt(n);
Fh = ifft((omega.' * ones(1, 2)) .* H) * sqrt(n);

dp = exp(1i * pi / n * (0 : 2 : 2*n-2));
dd = exp(1i * pi / n * (1 : 2 : 2*n));

C = hss('cauchy', dp, dd);
C = hss('low-rank', Gh, Fh) .* C;

fv = ifft(b);
fx = C \ fv;

x = omega' .* fft(fx);

end