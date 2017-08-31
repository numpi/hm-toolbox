function M = full(A)
n = A.ml+A.mr;
M = hss_mul(A, eye(n));
