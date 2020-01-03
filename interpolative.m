function [X, J] = interpolative(A, tol)

if isempty(A)
    J = []; X = A;
    return;
end

% compute interpolative (by columns) low rank approximation of A
if ~exist('tol', 'var')
    tol = hssoption('threshold');
end
[Q, R, p] = qr(A, 'vector');
for j = 1:length(p) % invert permutation
    ip(p(j))= j;
end
rk = sum(abs(diag(R))> tol * abs(R(1,1)));
J = p(1:rk);
X = R(1:rk, 1:rk)\ R(1:rk,ip);
end
