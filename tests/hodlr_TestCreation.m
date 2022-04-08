function hodlr_TestCreation
%
% Test the correction generation of hodlr matrices.
%

% Using a small block size hurts performances but helps to test out the
% code.
hodlroption('block-size', 32);
tol = hodlroption('threshold');

hnrm = @(A) norm(A, hodlroption('norm'));

switch hodlroption('norm')
    case 2
        C = @(A) log2(max(size(A)));
    case 'fro'
        C = @(A) sqrt(max(size(A)));
end

for n = [ 100, 1000 ]
	A = randn(n, n);

	H = hodlr(A);

	CheckTestResult(hnrm(A - full(H)), '<', C(A) * hnrm(A) * hodlroption('threshold'), ...
		'Generation of an hodlr representation for unstructured A');

	A = sprand(n, n, 4 / n);
	H = hodlr(A);

	CheckTestResult(hnrm(full(A) - full(H)), '<', C(A) * hnrm(full(A)) * hodlroption('threshold'), ...
		'Generation of an hodlr representation for a sparse A');
    
    A = zeros(n, n);
    H = hodlr(A);
    
    CheckTestResult(hnrm(A - full(H)), '<=', 0, ...
		'Generation of an hodlr representation for a zero matrix A');
end

for n = [ 100, 1000 ]
	x = linspace(0, 1, n);
	y = x + 1/(2*(n-1));
	A = 1 ./ (x - y.');

	H = hodlr(A);

	CheckTestResult(hnrm(A - full(H)), '<', 4 * C(A) * hnrm(A) * hodlroption('threshold'), ...
		'Generation of an hodlr representation for Cauchy A built from dense');
    
    H = hodlr('cauchy', -y, x);
    
    CheckTestResult(hnrm(A - full(H)), '<', 4 * C(A) * hnrm(A) * hodlroption('threshold'), ...
		'Generation of an hodlr representation for Cauchy A built with the Cauchy constructor');
end

n = 10000;
A = spdiags(ones(n, 1) * [-1 2 -1], -1:1, n, n);
H = hodlr('tridiagonal', A);

v = randn(n, 1);

CheckTestResult(norm(A*v - H*v), '<', normest(A, 1e-3) * norm(v) * tol, ...
    'Generation of an hodlr representation for tridiagonal A');

A = spdiags(randn(n, 1) * rand(1, 10), -4:5, n, n);
H = hodlr('banded', A);

v = randn(n, 1);

CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an hodlr representation for banded A');

k = 4;
U = rand(n, k); V = rand(n, k);

A = hodlr('low-rank', U, V);

CheckTestResult(norm(A*v - U*(V'*v)), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an hodlr representation for low-rank A');

n = 2048;

c = rand(10, 1); r = rand(1, 8); r(1) = c(1);
H = hodlr('toeplitz', c, r, n);

c(n) = 0; r(n) = 0;
A = toeplitz(c, r);

v = randn(n, 1);

CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an hodlr representation for banded Toeplitz A');

c = rand(n, 1) .* (.5.^(1 : n)'); r = rand(1, n) .* (.75.^(1 : n)); r(1) = c(1);
H = hodlr('toeplitz', c, r);

A = toeplitz(c, r);

CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an hodlr representation for Toeplitz A');

% Function samples
f = @(x,y) log(1 + abs(x - y));
x = linspace(0, 1, n);
H = hodlr('handle', @(i,j) f(x(j), x(i)'), n, n);

A = f( x, x' );

CheckTestResult(norm(A*v - H*v), '<', C(A) * hnrm(A) * norm(v) * tol, ...
    'Generation of an hodlr representation for A that samples f(x,y)');

% Generate random clustering
cluster = cumsum(randi(32, 1, 16));

n = cluster(end);

A = spdiags(randn(n, 1) * rand(1, 10), -4:5, n, n);
H = hodlr('banded', A, 'cluster', cluster);

v = randn(n, 1);

CheckTestResult(hnrm(A*v - H*v), '<', hnrm(v) * sqrt(n) * tol, ...
    'Generation of an hodlr representation for banded A (nonstandard clustering)');

k = 4;
U = rand(n, k); V = rand(n, k);

A = hodlr('low-rank', U, V, 'cluster', cluster);

CheckTestResult(hnrm(A*v - U*(V'*v)), '<', hnrm(v) * sqrt(n) * tol, ...
    'Generation of an hodlr representation for low-rank A (nonstandard clustering)');

c = rand(10, 1); r = rand(1, 8); r(1) = c(1);
H = hodlr('toeplitz', c, r, n, 'cluster', cluster);

c(n) = 0; r(n) = 0;
A = toeplitz(c, r);

v = randn(n, 1);

CheckTestResult(hnrm(A*v - H*v), '<', hnrm(v) * sqrt(n) * tol, ...
    'Generation of an hodlr representation for banded Toeplitz A (nonstandard clustering)');

c = rand(n, 1) .* (.5.^(1 : n)'); r = rand(1, n) .* (.75.^(1 : n)); r(1) = c(1);
H = hodlr('toeplitz', c, r, 'cluster', cluster);

A = toeplitz(c, r);

CheckTestResult(hnrm(A*v - H*v), '<', hnrm(v) * sqrt(n) * tol, ...
    'Generation of an hodlr representation for Toeplitz A (nonstandard clustering)');

lambda = rand(1,n); mu = rand(1,n);
V = randn(n, 4); R = randn(4, n);
L = randn(n, 3); W = randn(3, n);

H = hodlr('loewner', mu, lambda, V, R, L, W);
LO = (V*R - L*W) ./ (mu.' + lambda);

CheckTestResult(norm(full(H) - LO, 'fro'), '<', norm(LO, 'fro') * sqrt(n) * tol, ...
    'Generation of a hodlr representation for a Loewner matrix A');


H = hodlr('shifted-loewner', mu, lambda, V, R, L, W);
LO = (diag(mu)*V*R + L*W*diag(lambda)) ./ (mu.' + lambda);

CheckTestResult(norm(full(H) - LO, 'fro'), '<', norm(LO, 'fro') * sqrt(n) * tol, ...
    'Generation of a hodlr representation for a shifted Loewner matrix A');
end