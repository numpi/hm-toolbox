function hodlr_TestOperations
clear

hodlroption('block-size', 32);
hssoption('block-size', 32);
n = 512;

tol = hodlroption('threshold');
hnrm = @(A) norm(A, hodlroption('norm'));

A = randn(n);
H = hodlr(A)';
CheckTestResult(hnrm(A' - full(H)), '<', hnrm(A) * tol, ...
    'Transposition for unstructured hodlr');

A = tril(triu(randn(n), -4), 6);
H = hodlr('banded', A, 4, 6)';
CheckTestResult(hnrm(A' - full(H)), '<', hnrm(A) * tol, ...
    'Transposition for banded hodlr');

U = randn(n,7); V = randn(n,7);
H = hodlr('low-rank', U, V)';
CheckTestResult(hnrm(V * U' - full(H)), '<', hnrm(U) * hnrm(V) * tol, ...
    'Transposition for low-rank hodlr');

A = randn(n, n); hodlrA = hodlr(A); 
B = randn(n, n); hodlrB = hodlr(B); hssB = hss(B);

H = hodlr(A) + hodlr(B);
CheckTestResult(hnrm(A + B - full(H)), '<', hnrm(A) * tol, ...
    'Sum of the hodlr representation of unstructured A and unstructured B');

H = hodlrA + B;
CheckTestResult(hnrm(A + B - full(H)), '<', hnrm(A) * tol && isfloat(H), ...
    'Sum of the hodlr representation of  A and dense B');

H = A + hodlrB;
CheckTestResult(hnrm(A + B - full(H)), '<', hnrm(A) * tol && isfloat(H), ...
    'Sum of dense  A and hodlr representation of B');

H = hodlrA + hssB;
CheckTestResult(hnrm(A + B - full(H)), '<', hnrm(A) * tol && isa(H, 'hodlr'), ...
    'Sum of the hodlr representation of  A and HSS representation of  B');

H = hodlrA - B;
CheckTestResult(hnrm(A - B - full(H)), '<', hnrm(A) * tol && isfloat(H), ...
    'Difference of the hodlr representation of  A and dense B');

H = A - hodlrB;
CheckTestResult(hnrm(A - B - full(H)), '<', hnrm(A) * tol && isfloat(H), ...
    'Difference of dense  A and hodlr representation of B');

H = hodlrA - hssB;
CheckTestResult(hnrm(A - B - full(H)), '<', hnrm(A) * tol && isa(H, 'hodlr'), ...
    'Difference of the hodlr representation of  A and HSS representation of  B');


B = tril(triu(randn(n), -3),4);
H = hodlr(A) + hodlr('banded', B, 3, 4);
CheckTestResult(hnrm(A + B - full(H)), '<', hnrm(A) * tol, ...
    'Sum of the hodlr representation of unstructured A and banded B');

A = tril(triu(randn(n), -6),3);
H = hodlr('banded', A, 6, 3) + hodlr('banded', B, 3, 4);
CheckTestResult(hnrm(A + B - full(H)), '<', (hnrm(A) + hnrm(B)) * tol, ...
    'Sum of the hodlr representation of banded A and banded B');

U = randn(n, 6); V = randn(n, 6);
W = randn(n, 6); Z = randn(n, 6);
A = randn(n);
H = hodlr(A) + hodlr('low-rank', U, V);
CheckTestResult(hnrm(A + U*V' - full(H)), '<', (hnrm(A) + hnrm(U*V')) * tol, ...
    'Sum of the hodlr representation of unstructured A and low-rank B');

H = hodlr('low-rank', U, V) + hodlr('banded', B, 3, 4);
CheckTestResult(hnrm(U*V' + B - full(H)), '<', (hnrm(A) + hnrm(U*V')) * tol, ...
    'Sum of the hodlr representation of low-rank A and banded B');

H = hodlr('low-rank', U, V) + hodlr('low-rank', W, Z);
CheckTestResult(hnrm(U*V' + W*Z' - full(H)), '<', (hnrm(U*V') + hnrm(W*Z')) * tol, ...
    'Sum of the hodlr representation of low-rank A and low-rank B');

A = randn(n, n); hodlrA = hodlr(A); 
B = randn(n, n); hodlrB = hodlr(B); hssB = hss(B);
H = hodlr(A) * hodlr(B);
CheckTestResult(hnrm(A * B - full(H)), '<', hnrm(A) * hnrm(B) * tol, ...
    'Product of the hodlr representation of unstructured A and unstructured B');

H = hodlrA * B;
CheckTestResult(hnrm(A * B - full(H)), '<', hnrm(A) * tol && isfloat(H), ...
    'Product of the hodlr representation of  A and dense B');

H = A * hodlrB;
CheckTestResult(hnrm(A * B - full(H)), '<', hnrm(A) * tol && isfloat(H), ...
    'Product of dense  A and hodlr representation of B');

H = hodlrA * hssB;
CheckTestResult(hnrm(A * B - full(H)), '<', hnrm(A) * tol && isa(H, 'hodlr'), ...
    'Product of the hodlr representation of  A and HSS representation of  B');

B = tril(triu(randn(n), -3),4);
H = hodlr(A) * hodlr('banded', B, 3, 4);
CheckTestResult(hnrm(A * B - full(H)), '<', hnrm(A) * hnrm(B) * tol, ...
    'Product of the hodlr representation of unstructured A and banded B');

A = tril(triu(randn(n), -6),3);
H = hodlr('banded', A, 6, 3) * hodlr('banded', B, 3, 4);
CheckTestResult(hnrm(A * B - full(H)), '<', hnrm(A) * hnrm(B) * tol, ...
    'Product of the hodlr representation of banded A and banded B');

U = randn(n, 6); V = randn(n, 6);
W = randn(n, 6); Z = randn(n, 6);
A = randn(n);
H = hodlr(A) * hodlr('low-rank', U, V);
CheckTestResult(hnrm(A * U*V' - full(H)), '<', hnrm(A) * hnrm(U) * hnrm(V) * tol, ...
    'Product of the hodlr representation of unstructured A and low-rank B');

H = hodlr('low-rank', U, V) * hodlr('banded', B, 3, 4);
CheckTestResult(hnrm(U*V' * B - full(H)), '<', hnrm(U) * hnrm(V) * hnrm(B) * tol, ...
    'Product of the hodlr representation of low-rank A and banded B');

H = hodlr('low-rank', U, V) * hodlr('low-rank', W, Z);
CheckTestResult(hnrm(U*V' * W*Z' - full(H)), '<', hnrm(U) * hnrm(V) * hnrm(W) * hnrm(Z) * tol, ...
    'Product of the hodlr representation of low-rank A and low-rank B');

% Matrix vector multiplication
H = hodlrgallery('rand', n, 10);
x = rand(n, 3); b = H * x;
CheckTestResult(hnrm(b - full(H)*x), '<', hnrm(H) * hnrm(x) * eps, ...
    'Matrix vector multiplication');

% Hadamard product of random hodlr matrices
A = hodlrgallery('rand', n, 10);
B = hodlrgallery('rand', n, 10);
H = A .* B;
CheckTestResult(hnrm(full(A) .* full(B) - full(H)), '<', hnrm(full(A) .* full(B)) * tol, ...
    'Hadamard product of the hodlr representation of random hodlr A and  B');

A = randn(n, n); hodlrA = hodlr(A); 
B = randn(n, n); hodlrB = hodlr(B); hssB = hss(B);
H = hodlrA .* B;
CheckTestResult(hnrm(A .* B - full(H)), '<', hnrm(A) * tol && isfloat(H), ...
    'Hadamard product of the hodlr representation of  A and dense B');

H = A .* hodlrB;
CheckTestResult(hnrm(A .* B - full(H)), '<', hnrm(A) * tol && isfloat(H), ...
    'Hadamard product of dense  A and hodlr representation of B');

H = hodlrA .* hssB;
CheckTestResult(hnrm(A .* B - full(H)), '<', hnrm(A) * tol && isa(H, 'hodlr'), ...
    'Hadamard product of the hodlr representation of  A and HSS representation of  B');

A = randn(n, n); hodlrA = hodlr(A);
CheckTestResult(hnrm(tril(A) - tril(hodlrA)), '<', hnrm(A) * tol, ...
    'Tril of the hodlr representation of A ');

CheckTestResult(hnrm(triu(A) - full(triu(hodlrA))), '<', hnrm(A) * tol, ...
    'Triu of the hodlr representation of A ');

% Power of an hodlr matrix
A = hodlrgallery('rand', n, 10);
H = A^2;
CheckTestResult(hnrm(full(A)^2 - full(H)), '<', hnrm(full(A))^2 * tol, ...
    'Square of the hodlr representation of random hodlr A');

A = hodlrgallery('rand', n, 10);
H = A^3;
CheckTestResult(hnrm(full(A)^3 - full(H)), '<', hnrm(full(A))^3 * tol, ...
    'Cube of the hodlr representation of random hodlr A');

for p = [ 4, 7, 9, 12, 15 ]
    A = hodlrgallery('rand', n, 10);
    A = A / hnrm(A); % Make sure things stay bounded
    H = A^p;
    CheckTestResult(hnrm(full(A)^p - full(H)), '<', hnrm(full(A))^p * tol, ...
        sprintf('%d-th Power of the hodlr representation of random hodlr A', p));
end

% Linear systems
H = hodlrgallery('rand', n, 3);
x = rand(n, 5); b = H*x; y = H \ b;
CheckTestResult(hnrm(x - y), '<', cond(full(H)) * hnrm(b), ...
    'Solution of a linear system (implicit LU)');

[L, U] = lu(H); y = U \ (L \ b);
CheckTestResult(hnrm(x - y), '<', cond(full(H)) * hnrm(b), ...
    'Solution of a linear system (explicit LU)');



end
