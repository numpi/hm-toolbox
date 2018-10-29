function M = full(A)
n = size(A, 2);
if A.topnode == 0
    A.topnode = 1;
    M = hss_mul(A, eye(n));
    A.topnode = 0;
else
    M = hss_mul(A, eye(n));
end
