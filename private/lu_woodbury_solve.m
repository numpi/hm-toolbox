function y = lu_woodbury_solve(nu, mu, L, U, x, W, D, Z)
% As lu solve but the coefficient matrix has also the low-rank correction W * D * Z'
y = lu_solve(nu, mu, L, U, x);
t = lu_solve(nu, mu, L, U, x);
y = y - t * ((inv(D) + Z' * t) \ (Z' * y));
