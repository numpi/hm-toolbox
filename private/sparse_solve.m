function y = sparse_solve(nu, mu, L, U, p, q, x)

if nu > mu
    y = nu\(q*(U\(L\(p*x))));
else
    y = -mu \ x;
end

