function y = sparse_solve(nu, mu, L, U, p, q, x, A)

% [nu mu]
if nu > mu
    y = nu\(q*(U\(L\(p*x))));
else
    y = -mu \ x;
end


% [mu nu]

