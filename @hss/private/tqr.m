function [U, S, V] = tqr(A, tol)
%
% SVD-like interface to prrqr
%

[Q, R, p] = prrqr(A, tol, false, false);
% norm(Q*R - A(:,p))
U = Q;
ip(p) = 1 : length(p);
S = R(:, ip);

[V, S] = qr(S', 0); S = S';