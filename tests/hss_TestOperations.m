function hss_TestOperations
clear

hssoption('block-size', 32);
n = 512;

tol = hssoption('threshold');

A = randn(n);
H = hss(A)';
CheckTestResult(norm(A' - full(H)), '<', norm(A) * tol, ...
    'Transposition for unstructured HSS');

A = tril(triu(randn(n), -4), 6);
H = hss('banded', A, 4, 6)';
CheckTestResult(norm(A' - full(H)), '<', norm(A) * tol, ...
    'Transposition for banded HSS');

U = randn(n,7); V = randn(n,7);
H = hss('low-rank', U, V)';
CheckTestResult(norm(V * U' - full(H)), '<', norm(U) * norm(V) * tol, ...
    'Transposition for low-rank HSS');

A = randn(n, n);
B = randn(n, n);
H = hss(A) + hss(B);
CheckTestResult(norm(A + B - full(H)), '<', norm(A) * tol, ...
    'Sum of the HSS representation of unstructured A and unstructured B');

B = tril(triu(randn(n), -3),4);
H = hss(A) + hss('banded', B, 3, 4);
CheckTestResult(norm(A + B - full(H)), '<', norm(A) * tol, ...
    'Sum of the HSS representation of unstructured A and banded B');

A = tril(triu(randn(n), -6),3);
H = hss('banded', A, 6, 3) + hss('banded', B, 3, 4);
CheckTestResult(norm(A + B - full(H)), '<', norm(A) * tol, ...
    'Sum of the HSS representation of banded A and banded B');

U = randn(n, 6); V = randn(n, 6);
W = randn(n, 6); Z = randn(n, 6);
A = randn(n);
H = hss(A) + hss('low-rank', U, V);
CheckTestResult(norm(A + U*V' - full(H)), '<', norm(A) * tol, ...
    'Sum of the HSS representation of unstructured A and low-rank B');

H = hss('low-rank', U, V) + hss('banded', B, 3, 4);
CheckTestResult(norm(U*V' + B - full(H)), '<', norm(A) * tol, ...
    'Sum of the HSS representation of low-rank A and banded B');

H = hss('low-rank', U, V) + hss('low-rank', W, Z);
CheckTestResult(norm(U*V' + W*Z' - full(H)), '<', norm(A) * tol, ...
    'Sum of the HSS representation of low-rank A and low-rank B');

A = randn(n, n);
B = randn(n, n);
H = hss(A) * hss(B);
CheckTestResult(norm(A * B - full(H)), '<', norm(A) * norm(B) * tol, ...
    'Product of the HSS representation of unstructured A and unstructured B');

B = tril(triu(randn(n), -3),4);
H = hss(A) * hss('banded', B, 3, 4);
CheckTestResult(norm(A * B - full(H)), '<', norm(A) * norm(B) * tol, ...
    'Product of the HSS representation of unstructured A and banded B');

A = tril(triu(randn(n), -6),3);
H = hss('banded', A, 6, 3) * hss('banded', B, 3, 4);
CheckTestResult(norm(A * B - full(H)), '<', norm(A) * norm(B) * tol, ...
    'Product of the HSS representation of banded A and banded B');

U = randn(n, 6); V = randn(n, 6);
W = randn(n, 6); Z = randn(n, 6);
A = randn(n);
H = hss(A) * hss('low-rank', U, V);
CheckTestResult(norm(A * U*V' - full(H)), '<', norm(A) * norm(U) * norm(V) * tol, ...
    'Product of the HSS representation of unstructured A and low-rank B');

H = hss('low-rank', U, V) * hss('banded', B, 3, 4);
CheckTestResult(norm(U*V' * B - full(H)), '<', norm(U) * norm(V) * norm(B) * tol, ...
    'Product of the HSS representation of low-rank A and banded B');

H = hss('low-rank', U, V) * hss('low-rank', W, Z);
CheckTestResult(norm(U*V' * W*Z' - full(H)), '<', norm(U) * norm(V) * norm(W) * norm(Z) * tol, ...
    'Product of the HSS representation of low-rank A and low-rank B');

% Matrix vector multiplication
H = hssgallery('rand', n, 10);
x = rand(n, 3); b = H * x;
CheckTestResult(norm(b - full(H)*x), '<', norm(H) * norm(x) * eps, ...
    'Matrix vector multiplication');

% Hadamard product of random hss matrices
A = hssgallery('rand', n, 10);
B = hssgallery('rand', n, 10);
H = A .* B;
CheckTestResult(norm(full(A) .* full(B) - full(H)), '<', norm(full(A) .* full(B)) * tol, ...
    'Hadamard product of the HSS representation of random hss A and  B');
% Linear systems
H = hssgallery('rand', n, 3);
x = rand(n, 5); b = H*x; y = H \ b;
CheckTestResult(norm(x - y), '<', cond(full(H)) * norm(b), ...
    'Solution of a linear system (implicit ULV)');

F = ulv(H); y = ulv_solve(F, b);
CheckTestResult(norm(x - y), '<', cond(full(H)) * norm(b), ...
    'Solution of a linear system (explicit ULV)');



end
