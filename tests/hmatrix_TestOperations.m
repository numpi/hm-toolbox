function hmatrix_TestOperations
clear

hmatrixoption('block-size', 32);

n = 512;

tol = hmatrixoption('threshold');

A = randn(n);
H = hmatrix(A)';
CheckTestResult(norm(A' - full(H)), '<', norm(A) * tol, ...
    'Transposition for unstructured hmatrix');

A = tril(triu(randn(n), -4), 6);
H = hmatrix('banded', A, 4, 6)';
CheckTestResult(norm(A' - full(H)), '<', norm(A) * tol, ...
    'Transposition for banded hmatrix');

U = randn(n,7); V = randn(n,7);
H = hmatrix('low-rank', U, V)';
CheckTestResult(norm(V * U' - full(H)), '<', norm(U) * norm(V) * tol, ...
    'Transposition for low-rank hmatrix');

A = randn(n, n); hmatrixA = hmatrix(A);
B = randn(n, n); hmatrixB = hmatrix(B);

H = hmatrix(A) + hmatrix(B);
CheckTestResult(norm(A + B - full(H)), '<', norm(A) * tol, ...
    'Sum of the hmatrix representation of unstructured A and unstructured B');

B = tril(triu(randn(n), -3),4);
H = hmatrix(A) + hmatrix('banded', B, 3, 4);
CheckTestResult(norm(A + B - full(H)), '<', norm(A) * tol, ...
    'Sum of the hmatrix representation of unstructured A and banded B');

A = tril(triu(randn(n), -6),3);
H = hmatrix('banded', A, 6, 3) + hmatrix('banded', B, 3, 4);
CheckTestResult(norm(A + B - full(H)), '<', (norm(A) + norm(B)) * tol, ...
    'Sum of the hmatrix representation of banded A and banded B');

U = randn(n, 6); V = randn(n, 6);
W = randn(n, 6); Z = randn(n, 6);
A = randn(n);
H = hmatrix(A) + hmatrix('low-rank', U, V);
CheckTestResult(norm(A + U*V' - full(H)), '<', (norm(A) + norm(U*V')) * tol, ...
    'Sum of the hmatrix representation of unstructured A and low-rank B');

H = hmatrix('low-rank', U, V) + hmatrix('banded', B, 3, 4);
CheckTestResult(norm(U*V' + B - full(H)), '<', (norm(A) + norm(U*V')) * tol, ...
    'Sum of the hmatrix representation of low-rank A and banded B');

H = hmatrix('low-rank', U, V) + hmatrix('low-rank', W, Z);
CheckTestResult(norm(U*V' + W*Z' - full(H)), '<', (norm(U*V') + norm(W*Z')) * tol, ...
    'Sum of the hmatrix representation of low-rank A and low-rank B');

A = randn(n, n); hmatrixA = hmatrix(A);
B = randn(n, n); hmatrixB = hmatrix(B);
H = hmatrix(A) * hmatrix(B);
CheckTestResult(norm(A * B - full(H)), '<', norm(A) * norm(B) * tol, ...
    'Product of the hmatrix representation of unstructured A and unstructured B');




B = tril(triu(randn(n), -3),4);
H = hmatrix(A) * hmatrix('banded', B, 3, 4);
CheckTestResult(norm(A * B - full(H)), '<', norm(A) * norm(B) * tol, ...
    'Product of the hmatrix representation of unstructured A and banded B');

A = tril(triu(randn(n), -6),3);
H = hmatrix('banded', A, 6, 3) * hmatrix('banded', B, 3, 4);
CheckTestResult(norm(A * B - full(H)), '<', norm(A) * norm(B) * tol, ...
    'Product of the hmatrix representation of banded A and banded B');

U = randn(n, 6); V = randn(n, 6);
W = randn(n, 6); Z = randn(n, 6);
A = randn(n);
H = hmatrix(A) * hmatrix('low-rank', U, V);
CheckTestResult(norm(A * U*V' - full(H)), '<', norm(A) * norm(U) * norm(V) * tol, ...
    'Product of the hmatrix representation of unstructured A and low-rank B');

H = hmatrix('low-rank', U, V) * hmatrix('banded', B, 3, 4);
CheckTestResult(norm(U*V' * B - full(H)), '<', norm(U) * norm(V) * norm(B) * tol, ...
    'Product of the hmatrix representation of low-rank A and banded B');

H = hmatrix('low-rank', U, V) * hmatrix('low-rank', W, Z);
CheckTestResult(norm(U*V' * W*Z' - full(H)), '<', norm(U) * norm(V) * norm(W) * norm(Z) * tol, ...
    'Product of the hmatrix representation of low-rank A and low-rank B');

% Matrix vector multiplication
H = hmatrixgallery('rand', n, 10);
x = rand(n, 3); b = H * x;
CheckTestResult(norm(b - full(H)*x), '<', norm(H) * norm(x) * eps, ...
    'Matrix vector multiplication');

% Hadamard product of random hmatrix matrices
A = hmatrixgallery('rand', n, 10);
B = hmatrixgallery('rand', n, 10);
H = A .* B;
CheckTestResult(norm(full(A) .* full(B) - full(H)), '<', norm(full(A) .* full(B)) * tol, ...
    'Hadamard product of the hmatrix representation of random hmatrix A and  B');

A = randn(n, n); hmatrixA = hmatrix(A);
B = randn(n, n); hmatrixB = hmatrix(B);





A = randn(n, n); hmatrixA = hmatrix(A);
CheckTestResult(norm(tril(A) - full(tril(hmatrixA))), '<', norm(A) * tol, ...
    'Tril of the hmatrix representation of A ');

CheckTestResult(norm(triu(A) - full(triu(hmatrixA))), '<', norm(A) * tol, ...
    'Triu of the hmatrix representation of A ');

% Power of an hmatrix matrix
A = hmatrixgallery('rand', n, 10);
H = A^2;
CheckTestResult(norm(full(A)^2 - full(H)), '<', norm(full(A))^2 * tol, ...
    'Square of the hmatrix representation of random hmatrix A');

A = hmatrixgallery('rand', n, 10);
H = A^3;
CheckTestResult(norm(full(A)^3 - full(H)), '<', norm(full(A))^3 * tol, ...
    'Cube of the hmatrix representation of random hmatrix A');

for p = [ 4, 7, 9, 12, 15 ]
    A = hmatrixgallery('rand', n, 10);
    A = A / norm(A); % Make sure things stay bounded
    H = A^p;
    CheckTestResult(norm(full(A)^p - full(H)), '<', norm(full(A))^p * tol, ...
        sprintf('%d-th Power of the hmatrix representation of random hmatrix A', p));
end

% Linear systems
H = hmatrixgallery('rand', n, 3) + 2 * hmatrix('eye', n);
x = rand(n, 5); b = H*x; y = H \ b;
CheckTestResult(norm(x - y), '<', cond(full(H)) * norm(b), ...
    'Solution of a linear system (implicit LU)');

[L, U] = lu(H); y = U \ (L \ b);
CheckTestResult(norm(x - y), '<', cond(full(H)) * norm(b), ...
    'Solution of a linear system (explicit LU)');



end
