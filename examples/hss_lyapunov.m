%% Solving a Lyapunov equation in the HSS format.

%% Available solvers
% We consider the same problem of the <hm_lyapunov.html HODLR> case.
% However, the sign iteration is not available, since no inverse is
% implemented in the HSS format.
%
% On the other hand, an approach which often delivers better performance is
% presented in "Low-rank updates and a divide-and-conquer method for linear
% matrix equations", by D. Kressner, S. Massei, and L. Robol, 2018.
%
% This approach produces a solution in the HSS format, but the only
% ingredient needed is matrices with a fast mat-vec operation, and a fast
% linear system solution procedure. Matrices in the HSS format are a
% natural candidate, but sparse banded matrices work as well, and when the
% latter structure is available is usually more efficient.

n = 2^11;
h = 1 / (n + 1);

%% The solver in the HSS format
% We start by considering the pure HSS solver. We construct the same
% matrices that we had in the HODLR case.

x = linspace(0, 1, n + 2);
A = hss('banded', spdiags(ones(n, 1) * [ -1, 2, -1 ], -1:1, n, n), 1, 1);
C = hm2hss(hm('handle', @(i,j) log(1 + abs(x(j+1) - x(i+1)')), n, n));

%%
% Let us briefly check that all these matrices have a low-rank structure in
% their off-diagonal blocks:

fprintf('hssrank(A) = %d, hssrank(C) = %d\n', hssrank(A), hssrank(C));

%% Solving the equation
% We can now use LYAP to solve the Lyapunov equation. In this example, we
% use the sign iteration. This is also the default method used by the LYAP
% function.

tic; X = lyap(A, -h^2 * C); thss = toc;

%%
% Let us verify that the solution satifies the differential equation. Here,
% we compute the $L^2$ norm of the residual, relative to the $L^2$ norm of
% the solution $X$.
fprintf('[Standard D&C HSS solver] relative residual: %e, time: %fs\n', ...
    norm(A * X + X * A - h^2 * C, 'fro') / norm(X, 'fro'), ...
    thss);

%%
% Let's have a look at the solution
x = linspace(0, 1, n); mesh(x, x, full(X));

%% Using the sparse solver
% Since in this case we have a sparse representation of A, we can make use
% of it to speed up the computations.
sA = spdiags(ones(n, 1) * [-1 2 -1], -1:1, n, n);

tic; X = lyap(A, -h^2 * C, sA); tsparse = toc;

fprintf('[Sparse D&C HSS solver]   relative residual: %e, time: %fs\n', ...
    norm(A * X + X * A - h^2 * C, 'fro') / norm(X, 'fro'), ...
    tsparse);

%%
% Let's have a look at the solution
x = linspace(0, 1, n); mesh(x, x, full(X));



