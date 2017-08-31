function M = full(A)
n = size(A, 1);
M = hss_mul(A, eye(n));
