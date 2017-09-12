function hss_TestOperations
clear

	hssoption('block-size', 32);
	n = 512;

	A = randn(n);
	H = hss(A)';
	CheckTestResult(norm(A' - full(H)), '<', 1e3 * norm(A) * eps, ...
	'Transposition for unstructured HSS');

	A = tril(triu(randn(n), -4), 6);
	H = hss('banded', A, 4, 6)';
	CheckTestResult(norm(A' - full(H)), '<', 1e3 * norm(A) * eps, ...
	'Transposition for banded HSS');

	U = randn(n,7); V = randn(n,7);
	H = hss('low-rank', U, V)';
	CheckTestResult(norm(V * U' - full(H)), '<', 1e3 * norm(U) * norm(V) * eps, ...
	'Transposition for low-rank HSS');


	A = randn(n, n);
	B = randn(n, n);
	H = hss(A) + hss(B); 
	CheckTestResult(norm(A + B - full(H)), '<', 1e3 * norm(A) * eps, ...
	'Sum of the HSS representation of unstructured A and unstructured B');

	B = tril(triu(randn(n), -3),4);
	H = hss(A) + hss('banded', B, 3, 4); 
	CheckTestResult(norm(A + B - full(H)), '<', 1e3 * norm(A) * eps, ...
	'Sum of the HSS representation of unstructured A and banded B');	

	A = tril(triu(randn(n), -6),3);
	H = hss('banded', A, 6, 3) + hss('banded', B, 3, 4); 
	CheckTestResult(norm(A + B - full(H)), '<', 1e3 * norm(A) * eps, ...
	'Sum of the HSS representation of banded A and banded B');

	U = randn(n, 6); V = randn(n, 6);
	W = randn(n, 6); Z = randn(n, 6);
	A = randn(n);	
	H = hss(A) + hss('low-rank', U, V);
	CheckTestResult(norm(A + U*V' - full(H)), '<', 1e3 * norm(A) * eps, ...
	'Sum of the HSS representation of unstructured A and low-rank B');

	H = hss('low-rank', U, V) + hss('banded', B, 3, 4); 
	CheckTestResult(norm(U*V' + B - full(H)), '<', 1e3 * norm(A) * eps, ...
	'Sum of the HSS representation of low-rank A and banded B');

	H = hss('low-rank', U, V) + hss('low-rank', W, Z); 
	CheckTestResult(norm(U*V' + W*Z' - full(H)), '<', 1e3 * norm(A) * eps, ...
	'Sum of the HSS representation of low-rank A and low-rank B');

	A = randn(n, n);
	B = randn(n, n);
	H = hss(A) * hss(B); 
	CheckTestResult(norm(A * B - full(H)), '<', 1e3 * norm(A) * norm(B) * eps, ...
	'Product of the HSS representation of unstructured A and unstructured B');

	B = tril(triu(randn(n), -3),4);
	H = hss(A) * hss('banded', B, 3, 4); 
	CheckTestResult(norm(A * B - full(H)), '<', 1e3 * norm(A) * norm(B) * eps, ...
	'Product of the HSS representation of unstructured A and banded B');	

	A = tril(triu(randn(n), -6),3);
	H = hss('banded', A, 6, 3) * hss('banded', B, 3, 4); 
	CheckTestResult(norm(A * B - full(H)), '<', 1e3 * norm(A) * norm(B) * eps, ...
	'Product of the HSS representation of banded A and banded B');

	U = randn(n, 6); V = randn(n, 6);
	W = randn(n, 6); Z = randn(n, 6);
	A = randn(n);	
	H = hss(A) * hss('low-rank', U, V);
	CheckTestResult(norm(A * U*V' - full(H)), '<', 1e3 * norm(A) * norm(U) * norm(V) * eps, ...
	'Product of the HSS representation of unstructured A and low-rank B');

	H = hss('low-rank', U, V) * hss('banded', B, 3, 4); 
	CheckTestResult(norm(U*V' * B - full(H)), '<', 1e3 * norm(U) * norm(V) * norm(B) * eps, ...
	'Product of the HSS representation of low-rank A and banded B');

	H = hss('low-rank', U, V) * hss('low-rank', W, Z); 
	CheckTestResult(norm(U*V' * W*Z' - full(H)), '<', 1e3 * norm(U) * norm(V) * norm(W) * norm(Z) * eps, ...
	'Product of the HSS representation of low-rank A and low-rank B');

end
