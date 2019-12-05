function y = sparse_solve(nu, mu, L, U, p, q, x)

if nu > mu
    y = nu\((U\(L\(x(p,:))))); y = y(q, :);
else
    y = -mu \ x;
end

