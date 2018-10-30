function M = full(A)

m = size(A, 1);
n = size(A, 2);

if m < n
    M = full(A')';
    return;
end

if A.topnode == 0
    A.topnode = 1;
    M = hss_mul(A, eye(n));
    A.topnode = 0;
else
    M = hss_mul(A, eye(n));
end
