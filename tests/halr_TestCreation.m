function halr_TestCreation
%
% Test the correction generation of halr matrices.
%

% Using a small block size hurts performances but helps to test out the
% code.
halroption('block-size', 32);
tol = halroption('threshold');

for n = [ 100, 1000 ]
    A = randn(n, n);
    
    H = halr(A);
    
    CheckTestResult(norm(A - full(H)), '<', norm(A) * halroption('threshold'), ...
        'Generation of an halr representation for unstructured A');
    
    A = sprand(n, n, 4 / n);
    H = halr(A);
    
    CheckTestResult(norm(full(A) - full(H)), '<', norm(full(A)) * halroption('threshold'), ...
        'Generation of an halr representation for a sparse A');
end

for n = [ 100, 1000 ]
    x = linspace(0, 1, n);
    y = x + 1/(2*(n-1));
    A = 1 ./ (x - y.');
    
    H = halr(A);
    
    CheckTestResult(norm(A - full(H)), '<', 10 * norm(A) * halroption('threshold'), ...
        'Generation of an halr representation for Cauchy A built from dense');
    
    %H = halr('cauchy', -y, x);
    
    %CheckTestResult(norm(A - full(H)), '<', 10 * log2(n) * norm(A) * halroption('threshold'), ...
    %	'Generation of an halr representation for Cauchy A built with the Cauchy constructor');
end

n = 10000;
% A = spdiags(ones(n, 1) * [-1 2 -1], -1:1, n, n);
% H = halr('banded', A);
%
% v = randn(n, 1);
%
% CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * eps, ...
%     'Generation of an halr representation for tridiagonal A');
%
% A = spdiags(randn(n, 1) * rand(1, 10), -4:5, n, n);
% H = halr('banded', A);
%
% v = randn(n, 1);
%
% CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * tol, ...
%     'Generation of an halr representation for banded A');
%
% k = 4;
% U = rand(n, k); V = rand(n, k);
%
% A = halr('low-rank', U, V);
%
% CheckTestResult(norm(A*v - U*(V'*v)), '<', norm(v) * sqrt(n) * tol, ...
%     'Generation of an halr representation for low-rank A');
%
n = 2048;
%
% c = rand(10, 1); r = rand(1, 8); r(1) = c(1);
% H = halr('toeplitz', c, r, n);
%
% c(n) = 0; r(n) = 0;
% A = toeplitz(c, r);
%
v = randn(n, 1);
%
% CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * tol, ...
%     'Generation of an halr representation for banded Toeplitz A');
%
% c = rand(n, 1) .* (.5.^(1 : n)'); r = rand(1, n) .* (.75.^(1 : n)); r(1) = c(1);
% H = halr('toeplitz', c, r);
%
% A = toeplitz(c, r);
%
% CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * tol, ...
%     'Generation of an halr representation for Toeplitz A');

% Function samples
f = @(x,y) log(1 + abs(x - y));
x = linspace(0, 1, n);
H = halr('handle', @(i,j) f(x(j), x(i)'), n, n);

A = f( x, x' );

CheckTestResult(norm(A*v - H*v), '<', norm(A) * norm(v) * sqrt(n) * halroption('threshold'), ...
    'Generation of an halr representation for A that samples f(x,y)');

return

% Generate random clustering
cluster = cumsum(randi(32, 1, 16));
n = cluster(end);

A = spdiags(randn(n, 1) * rand(1, 10), -4:5, n, n);
H = halr('banded', A, 'cluster', cluster);

v = randn(n, 1);

CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an halr representation for banded A (nonstandard clustering)');

k = 4;
U = rand(n, k); V = rand(n, k);

A = halr('low-rank', U, V, 'cluster', cluster);

CheckTestResult(norm(A*v - U*(V'*v)), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an halr representation for low-rank A (nonstandard clustering)');

c = rand(10, 1); r = rand(1, 8); r(1) = c(1);
H = halr('toeplitz', c, r, n, 'cluster', cluster);

c(n) = 0; r(n) = 0;
A = toeplitz(c, r);

v = randn(n, 1);

CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an halr representation for banded Toeplitz A (nonstandard clustering)');

c = rand(n, 1) .* (.5.^(1 : n)'); r = rand(1, n) .* (.75.^(1 : n)); r(1) = c(1);
H = halr('toeplitz', c, r, 'cluster', cluster);

A = toeplitz(c, r);

CheckTestResult(norm(A*v - H*v), '<', norm(v) * sqrt(n) * tol, ...
    'Generation of an halr representation for Toeplitz A (nonstandard clustering)');


end
