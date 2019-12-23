function hss_TestCreation
%
% Test the correction generation of HSS matrices.
%

% Using a small block size hurts performances but helps to test out the
% code.
clear
hssoption('block-size', 32);
hnrm = @(A) norm(A, hssoption('norm'));
C = @(A) sqrt(max(size(A)));
    

tol = hssoption('threshold');

n = 512;
A = randn(n, n);

H = hss(A);

CheckTestResult(norm(A - full(H)), '<', C(A) * hnrm(A) * tol, ...
    'Generation of an HSS representation for unstructured A');

n = 2^15;
A = spdiags(randn(n, 6) , -2:3, n, n);
H = hss('banded', A, 2, 3);

v = randn(n, 1);

CheckTestResult(norm(A*v - H*v), '<', C(A) * hnrm(H) * tol * norm(v), ...
    'Generation of an HSS representation for banded A');

n = 2^15;
U = randn(n,5);
V = randn(n,5);
B = randn(5);
v = randn(n, 10);

H = hss('low-rank', U, V);

CheckTestResult(norm(U * (V'*v) - H*v), '<', C(A) * hnrm(H) * tol * norm(v), ...
    'Generation of an HSS representation for U * V^T ');

n = 2^15;
U = randn(n,5);
V = randn(n,5);
B = randn(5);
v = randn(n, 10);

H = hss('low-rank', U, V, B);

CheckTestResult(norm(U * B *(V'*v) - H*v), '<', C(A) * hnrm(H) * tol * norm(v), ...
    'Generation of an HSS representation for U * B * V^T ');

n = 1024;
k = 3;
H = hss(diag(rand(n, 1)) + tril(rand(n, k) * rand(k, n)) + triu(rand(n, k) * rand(k, n)));

CheckTestResult(hssrank(H), '<=', 2*k, ...
    'Checking the rank of dense matrices (dense constructor)');

for n = [ 100, 1000 ]
	x = linspace(0, 1, n);
	y = x + 1/(2*(n-1));
	A = 1 ./ (x - y.');

	H = hss(A);

	CheckTestResult(norm(A - full(H)), '<', C(A) * hnrm(A) * tol, ...
		'Generation of an HSS representation for Cauchy A built from dense');
    
    H = hss('cauchy', -y, x);
    
    CheckTestResult(norm(A - full(H)), '<', C(A) * hnrm(A) * tol, ...
		'Generation of an HSS representation for Cauchy A built with the Cauchy constructor');
end


end
