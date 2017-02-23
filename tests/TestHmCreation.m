function TestHmCreation
%
% Test the correction generation of HM matrices. 
%

% Using a small block size hurts performances but helps to test out the
% code. 
hmoption('block-size', 32);

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
	'Generation of an HM representation for banded A');

end