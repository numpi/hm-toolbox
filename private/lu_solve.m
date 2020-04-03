function y = lu_solve(nu, mu, L, U, x)

if nu > mu
    y = nu\(U\(L\x));
else
    y = -mu \ x;
end

