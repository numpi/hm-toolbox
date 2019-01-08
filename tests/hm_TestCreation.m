function hm_TestCreation
%
% Test the correction generation of HM matrices.
%

% Using a small block size hurts performances but helps to test out the
% code.
hmoption('block-size', 32);
tol = hmoption('threshold');

n = 100;
A = randn(n, n);

H = hm(A);

CheckTestResult(norm(A - full(H)), '<', 1e3 * norm(A) * eps, ...
    'Generation of an HM representation for unstructured A');

n = 10000;
A = spdiags(ones(n, 1) * [-1 2 -1], -1:1, n, n);
H = hm('tridiagonal', A);

v = randn(n, 1);

CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * eps, ...
    'Generation of an HM representation for tridiagonal A');

A = spdiags(randn(n, 1) * rand(1, 10), -4:5, n, n);
H = hm('banded', A);

v = randn(n, 1);

CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an HM representation for banded A');

k = 4;
U = rand(n, k); V = rand(n, k);

A = hm('low-rank', U, V);

CheckTestResult(norm(A*v - U*(V'*v)), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an HM representation for low-rank A');

n = 2048;

c = rand(10, 1); r = rand(1, 8); r(1) = c(1);
H = hm('toeplitz', c, r, n);

c(n) = 0; r(n) = 0;
A = toeplitz(c, r);

v = randn(n, 1);

CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an HM representation for banded Toeplitz A');

c = rand(n, 1) .* (.5.^(1 : n)'); r = rand(1, n) .* (.75.^(1 : n)); r(1) = c(1);
H = hm('toeplitz', c, r);

A = toeplitz(c, r);

CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an HM representation for Toeplitz A');

% Chebfun2
f = @(x,y) log(1 + abs(x - y));
H = hm('chebfun2', f, [0 1], [0 1], n, n);

A = f( linspace(0, 1, n), linspace(0, 1, n)' );

CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an HM representation for A that samples f(x,y)');

% Generate random clustering
cluster = cumsum(randi(32, 1, 16));

n = cluster(end);

A = spdiags(randn(n, 1) * rand(1, 10), -4:5, n, n);
H = hm('banded', A, 'cluster', cluster);

v = randn(n, 1);

CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an HM representation for banded A (nonstandard clustering)');

k = 4;
U = rand(n, k); V = rand(n, k);

A = hm('low-rank', U, V, 'cluster', cluster);

CheckTestResult(norm(A*v - U*(V'*v)), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an HM representation for low-rank A (nonstandard clustering)');

c = rand(10, 1); r = rand(1, 8); r(1) = c(1);
H = hm('toeplitz', c, r, n, 'cluster', cluster);

c(n) = 0; r(n) = 0;
A = toeplitz(c, r);

v = randn(n, 1);

CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an HM representation for banded Toeplitz A (nonstandard clustering)');

c = rand(n, 1) .* (.5.^(1 : n)'); r = rand(1, n) .* (.75.^(1 : n)); r(1) = c(1);
H = hm('toeplitz', c, r, 'cluster', cluster);

A = toeplitz(c, r);

CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an HM representation for Toeplitz A (nonstandard clustering)');


end