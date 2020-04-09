function y = sparse_woodbury_solve(nu, mu, L, U, p, q, x, W, Z)
% As sparse solve but the coefficient matrix has also the low-rank correction W * Z'
y = sparse_solve(nu, mu, L, U, p, q, x); 
t = sparse_solve(nu, mu, L, U, p, q, W);
y = y - t * ((eye(size(W, 2)) + Z' * t) \ (Z' * y)); 

