%% Solving a rank structured linear system
%
% Solve a linear system defined by the sum of a banded matrix with a
% low-rank one, using the HODLR representation.

n = 16384;

fprintf('\nHM_LINEAR_SYSTEM // Solve a linear system with %d unknowns\n\n', n);

% We construct a banded matrix, and make sure it is invertible by shifting
% the diagonal entries.
spA = spdiags(rand(n, 3) + 2 * ones(n, 1) * [ 0 1 0 ], -1 : 1, n, n);
bandA = 1; % Bandwidth of A

% We construct its HODLR representation, and sum a random low-rank matrix
U = rand(n, 2); V = rand(n, 2);
A = hm('banded', spA, bandA) + hm('low-rank', U, V);

% Construct a random RHS ans solve the linear system
b = rand(n, 1);
x = A \ b;

r = norm(A * x - b);
fprintf(' - Residue of the linear system: %e\n', r);

if exist('timeit')
    fprintf('   Time needed to solve the system: %1.3f seconds\n\n', ...
        timeit(@() A \ b, 1));
end

% If more solutions are needed, it might be useful to precompute the LU
% factorization.
[L, U] = lu(A);
x = U \ (L \ b);

if exist('timeit')
    fprintf(' - Time needed to solve the system by backsubstution: %1.3f seconds\n', ...
        timeit(@() U \ (L \ b), 1));
    fprintf('   Time needed for the LU: %1.3f seconds\n', timeit(@() lu(A), 2));
end

