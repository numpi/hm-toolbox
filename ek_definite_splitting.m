function [Up, Un] = ek_definite_splitting(U, V, tol)
%DEFINITE_SPLITTING

[Q, R] = qr(U, 0);

W = R * (V' * Q);

% Ensure that W is symmetric
W = .5 * (W + W');

% Diagonalize W
[Z,D] = eig(W);

d = diag(D);

nrm = max(abs(d));

pos_indices = find(d > tol * nrm);
neg_indices = find(d < -tol * nrm);

Up = Q * Z(:,pos_indices) * diag(sqrt(d(pos_indices)));
Un = Q * Z(:,neg_indices) * diag(sqrt(-d(neg_indices)));

end

