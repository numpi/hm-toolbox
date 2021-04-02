function halr_TestOperations
clear

halroption('block-size', 32);

n = 512;

tol = halroption('threshold');

A = randn(n);
H = halr(A)';
CheckTestResult(norm(A' - full(H)), '<', norm(A) * tol, ...
    'Transposition for unstructured halr');

A = tril(triu(randn(n), -4), 6);
H = halr('banded', A, 4, 6)';
CheckTestResult(norm(A' - full(H)), '<', norm(A) * tol, ...
    'Transposition for banded halr');

U = randn(n,7); V = randn(n,7);
H = halr('low-rank', U, V)';
CheckTestResult(norm(V * U' - full(H)), '<', norm(U) * norm(V) * tol, ...
    'Transposition for low-rank halr');

A = randn(n, n); halrA = halr(A);
B = randn(n, n); halrB = halr(B);

H = halr(A) + halr(B);
CheckTestResult(norm(A + B - full(H)), '<', norm(A) * tol, ...
    'Sum of the halr representation of unstructured A and unstructured B');

B = tril(triu(randn(n), -3),4);
H = halr(A) + halr('banded', B, 3, 4);
CheckTestResult(norm(A + B - full(H)), '<', norm(A) * tol, ...
    'Sum of the halr representation of unstructured A and banded B');

A = tril(triu(randn(n), -6),3);
H = halr('banded', A, 6, 3) + halr('banded', B, 3, 4);
CheckTestResult(norm(A + B - full(H)), '<', (norm(A) + norm(B)) * tol, ...
    'Sum of the halr representation of banded A and banded B');

U = randn(n, 6); V = randn(n, 6);
W = randn(n, 6); Z = randn(n, 6);
A = randn(n);
H = halr(A) + halr('low-rank', U, V);
CheckTestResult(norm(A + U*V' - full(H)), '<', (norm(A) + norm(U*V')) * tol, ...
    'Sum of the halr representation of unstructured A and low-rank B');

H = halr('low-rank', U, V) + halr('banded', B, 3, 4);
CheckTestResult(norm(U*V' + B - full(H)), '<', (norm(A) + norm(U*V')) * tol, ...
    'Sum of the halr representation of low-rank A and banded B');

H = halr('low-rank', U, V) + halr('low-rank', W, Z);
CheckTestResult(norm(U*V' + W*Z' - full(H)), '<', (norm(U*V') + norm(W*Z')) * tol, ...
    'Sum of the halr representation of low-rank A and low-rank B');

A = randn(n, n); halrA = halr(A);
B = randn(n, n); halrB = halr(B);
H = halrA * halrB;
CheckTestResult(norm(A * B - full(H)), '<', norm(A) * norm(B) * tol, ...
    'Product of the halr representation of unstructured A and unstructured B');

B = tril(triu(randn(n), -3),4);
H = halr(A) * halr('banded', B, 3, 4);
CheckTestResult(norm(A * B - full(H)), '<', norm(A) * norm(B) * tol, ...
    'Product of the halr representation of unstructured A and banded B');

A = tril(triu(randn(n), -6),3);
H = halr('banded', A, 6, 3) * halr('banded', B, 3, 4);
CheckTestResult(norm(A * B - full(H)), '<', norm(A) * norm(B) * tol, ...
    'Product of the halr representation of banded A and banded B');

U = randn(n, 6); V = randn(n, 6);
W = randn(n, 6); Z = randn(n, 6);
A = randn(n);
H = halr(A) * halr('low-rank', U, V);
CheckTestResult(norm(A * U*V' - full(H)), '<', norm(A) * norm(U) * norm(V) * tol, ...
    'Product of the halr representation of unstructured A and low-rank B');

H = halr('low-rank', U, V) * halr('banded', B, 3, 4);
CheckTestResult(norm(U*V' * B - full(H)), '<', norm(U) * norm(V) * norm(B) * tol, ...
    'Product of the halr representation of low-rank A and banded B');

H = halr('low-rank', U, V) * halr('low-rank', W, Z);
CheckTestResult(norm(U*V' * W*Z' - full(H)), '<', norm(U) * norm(V) * norm(W) * norm(Z) * tol, ...
    'Product of the halr representation of low-rank A and low-rank B');

% Matrix vector multiplication
H = halrgallery('rand', n, 10);
x = rand(n, 3); b = H * x;
CheckTestResult(norm(b - full(H)*x), '<', norm(H) * norm(x) * eps, ...
    'Matrix vector multiplication');

% Hadamard product of random halr matrices
A = halrgallery('rand', n, 10);
B = halrgallery('rand', n, 10);
H = A .* B;
CheckTestResult(norm(full(A) .* full(B) - full(H)), '<', norm(full(A) .* full(B)) * tol, ...
    'Hadamard product of the halr representation of random halr A and  B');

A = randn(n, n); halrA = halr(A);
CheckTestResult(norm(tril(A) - full(tril(halrA))), '<', norm(A) * tol, ...
    'Tril of the halr representation of A ');

CheckTestResult(norm(triu(A) - full(triu(halrA))), '<', norm(A) * tol, ...
    'Triu of the halr representation of A ');

% Power of an halr matrix
A = halrgallery('rand', n, 10);
H = A^2;
CheckTestResult(norm(full(A)^2 - full(H)), '<', norm(full(A))^2 * tol, ...
    'Square of the halr representation of random halr A');

A = halrgallery('rand', n, 10);
H = A^3;
CheckTestResult(norm(full(A)^3 - full(H)), '<', norm(full(A))^3 * tol, ...
    'Cube of the halr representation of random halr A');

for p = [ 4, 7, 9, 12, 15 ]
    A = halrgallery('rand', n, 10);
    A = A / norm(A); % Make sure things stay bounded
    H = A^p;
    CheckTestResult(norm(full(A)^p - full(H)), '<', norm(full(A))^p * tol, ...
        sprintf('%d-th Power of the halr representation of random halr A', p));
end

% Linear systems
H = halrgallery('rand', n, 3) + 2 * halr('eye', n);
x = rand(n, 5); b = H*x; y = H \ b;
CheckTestResult(norm(x - y), '<', cond(full(H)) * norm(b), ...
    'Solution of a linear system (implicit LU)');

[L, U] = lu(H); y = U \ (L \ b);
CheckTestResult(norm(x - y), '<', cond(full(H)) * norm(b), ...
    'Solution of a linear system (explicit LU)');



end
