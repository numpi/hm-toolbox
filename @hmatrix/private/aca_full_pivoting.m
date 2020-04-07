function [U, V, I, J] = aca_full_pivoting(Afun, m, n, tol, comp, debug)
%
% Construct a low-rank representation of Afun
% by means of adaptive cross approximation with full pivoting (heuristic)
%
% ----------------------INPUT--------------------------------------------------
%
% Afun   matrix or function handle returning entries of the matrix
% m, n	 rows and cols of the matrix to compress
% tol    stopping criterion for the Adaptive Cross Approximation
% comp   re-compression at the end (optional), default = 0
% debug  enables some debugging prints (optional), default = 0
%
% ----------------------Output-------------------------------------------------
%
% U,V    factors of the low-rank approximation Afun = U * V'
% I	 row indices of the cross
% J	 column indices of the cross
% -----------------------------------------------------------------------------

if ~exist('comp', 'var')
    comp = 0;
end
if ~exist('debug', 'var')
    debug = 0;
end

U = zeros(m, 0);
V = zeros(n, 0);
I = [];
J = [];
k = 1;

while k <= min(m,n)
    % Compute the residual matrix
    R = Afun(1:m, 1:n) - U * V';
    % Select the pivot as maximal element in the residual matrix
    [temp, indj] = max(abs(R), [], 2);
    [~, indi] = max(temp);
    indj = indj(indi);
    piv = R(indi, indj);

    % Compute new factors
    a = R((1:m)', indj);
    b = R(indi, 1:n);
    a = a/piv;
    U = [U, a];  I = [I, indi];
    V = [V, b']; J = [J, indj];

    % Stopping criterion (heuristic)
    if k == 1
        nrm = norm(a) * norm(b);
    end
    if norm(a) * norm(b) < tol * nrm
        break
    end
    if debug
        fprintf('Iteration: %d Pivot at (%d,%d): %e\n', k, ind, new_ind, b(new_ind))
    end
    k = k + 1;
end
if comp
    [QU, RU] = qr(U, 0); % Re-compression
    [QV, RV] = qr(V, 0);
    [U, S, V] = svd(RU * RV', 'econ');
    rk = sum(diag(S) > tol * S(1,1));
    U = QU * U(:,1:rk) * sqrt(S(1:rk,1:rk));
    V = QV * V(:,1:rk) * sqrt(S(1:rk,1:rk));
end

