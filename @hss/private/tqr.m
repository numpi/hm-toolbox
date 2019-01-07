function [U, S, V] = tqr(A,tol)
if min(size(A)) == 0
    U = zeros(size(A, 1), 0); S = []; V = zeros(size(A, 2), 0); return;
end
[U, R, p] = qr(A, 0);
ip = zeros(length(p),1);
for i = 1:length(ip)
    ip(p(i)) = i;
end
n = sum(diag(abs(R)) > tol);
U = U(:, 1:n);
[S, V] = qr(R(1:n, ip), 0);
V = V';


