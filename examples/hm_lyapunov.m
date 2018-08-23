%% Solving a Lyapunov equation

%% Problem setting
%
% We consider the problem $AX + XA^T = C$, with $A$ symmetric positive
% definite, and $C$ with a HODLR structure. This can be solved by using the
% LYAP function in the HM format. If you are interested in the solver in
% the HSS format, check this <hss_lyapunov.html page>.
%
% Assume that $A = {\mathrm{diag}}_n(-1, 2, -1) \in {C}^{n \times n}$,
% and $C$ is the matrix with the sampling of $g(x,y) = \log(1 + |x - y|)$,
% scaled by a factor $h^2$, with $h = \frac{1}{n + 1}$.
%
% Then, the solution $X$ is a sampling of the function $u(x,y)$ that solves
% the PDE with zero Dirichlet boundary conditions
%
% $$
%   \frac{\partial^2}{\partial x^2} u(x,y) +
%    \frac{\partial^2}{\partial y^2} u(x,y) = \log(1 + |x - y|),
% $$
%
% where $(x,y) \in [0, 1]^2$, and $u(x,y) \equiv 0$ on $\partial [0, 1]^2$.
%
% For this test, we select $n = 2^{11}$.

n = 2^11;
h = 1 / (n + 1);

%% Building the matrices
% We start by constructing the matrices. For $A$, we can use the banded
% constructor:

A = hm('banded', spdiags(ones(n, 1) * [ -1, 2, -1 ], -1:1, n, n), 1);

% To build C, instead, we rely on chebfun2 (make sure to have chebfun in
% your MATLAB path for the following command to work!).

% C = hm('chebfun2', @(x,y) log(1 + abs(x - y)), [0 1], [0 1], n);
C = hm('chebfun2', @(x,y) log(1 + abs(x - y)), [0 1], [0 1], n);

% Let us briefly check that all these matrices have a low-rank structure in
% their off-diagonal blocks:

[ hmrank(A), hmrank(C) ]

%% Solving the equation
% We can now use LYAP to solve the Lyapunov equation. In this example, we
% use the sign iteration. This is also the default method used by the LYAP
% function.

tic; X = lyap(A, -h^2 * C, 'method', 'sign'); toc

%%
% Let us verify that the solution satifies the differential equation. Here,
% we compute the $L^2$ norm of the residual, relative to the $L^2$ norm of
% the solution $X$.
norm(A * X + X * A - h^2 * C, 'fro') / norm(X, 'fro')

%%
% Let's have a look at the solution
x = linspace(0, 1, n); mesh(x, x, full(X));

%% Other solution methods
% Alternative solvers are available. For instance, we have one based on the
% integral formula for the solution:
%
% $$
%   X = \int_0^\infty e^{-tA} C e^{-tA} dt.
% $$.
%
% Most of the time, the solver based on the sign iteration will outperform
% this option in speed and accuracy. However, this approach is easily
% parallelizable, and can this can be leveraged by adding the parameters
% 'parallel', true to the function call. Note that if you see these pages
% online, they probably have been generated where no parallelism was
% available, so expect this approach to be slower.

tic; X = lyap(A, -h^2 * C, 'method', 'expm', 'parallel', true); toc
norm(A * X + X * A - h^2 * C, 'fro') / norm(X, 'fro')

%%
% Let's have a look at the solution
x = linspace(0, 1, n); mesh(x, x, full(X));
