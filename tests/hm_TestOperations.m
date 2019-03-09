function hm_TestOperations
clear

hmoption('block-size', 32);
hssoption('block-size', 32);
n = 512;

tol = hmoption('threshold');

A = randn(n);
H = hm(A)';
CheckTestResult(norm(A' - full(H)), '<', norm(A) * tol, ...
    'Transposition for unstructured hm');

A = tril(triu(randn(n), -4), 6);
H = hm('banded', A, 4, 6)';
CheckTestResult(norm(A' - full(H)), '<', norm(A) * tol, ...
    'Transposition for banded hm');

U = randn(n,7); V = randn(n,7);
H = hm('low-rank', U, V)';
CheckTestResult(norm(V * U' - full(H)), '<', norm(U) * norm(V) * tol, ...
    'Transposition for low-rank hm');

A = randn(n, n); hmA = hm(A); 
B = randn(n, n); hmB = hm(B); hssB = hss(B);

H = hm(A) + hm(B);
CheckTestResult(norm(A + B - full(H)), '<', norm(A) * tol, ...
    'Sum of the HM representation of unstructured A and unstructured B');

H = hmA + B;
CheckTestResult(norm(A + B - full(H)), '<', norm(A) * tol && isfloat(H), ...
    'Sum of the HM representation of  A and dense B');

H = A + hmB;
CheckTestResult(norm(A + B - full(H)), '<', norm(A) * tol && isfloat(H), ...
    'Sum of dense  A and HM representation of B');

H = hmA + hssB;
CheckTestResult(norm(A + B - full(H)), '<', norm(A) * tol && isa(H, 'hm'), ...
    'Sum of the HM representation of  A and HSS representation of  B');

H = hmA - B;
CheckTestResult(norm(A - B - full(H)), '<', norm(A) * tol && isfloat(H), ...
    'Difference of the HM representation of  A and dense B');

H = A - hmB;
CheckTestResult(norm(A - B - full(H)), '<', norm(A) * tol && isfloat(H), ...
    'Difference of dense  A and HM representation of B');

H = hmA - hssB;
CheckTestResult(norm(A - B - full(H)), '<', norm(A) * tol && isa(H, 'hm'), ...
    'Difference of the HM representation of  A and HSS representation of  B');


B = tril(triu(randn(n), -3),4);
H = hm(A) + hm('banded', B, 3, 4);
CheckTestResult(norm(A + B - full(H)), '<', norm(A) * tol, ...
    'Sum of the hm representation of unstructured A and banded B');

A = tril(triu(randn(n), -6),3);
H = hm('banded', A, 6, 3) + hm('banded', B, 3, 4);
CheckTestResult(norm(A + B - full(H)), '<', (norm(A) + norm(B)) * tol, ...
    'Sum of the hm representation of banded A and banded B');

U = randn(n, 6); V = randn(n, 6);
W = randn(n, 6); Z = randn(n, 6);
A = randn(n);
H = hm(A) + hm('low-rank', U, V);
CheckTestResult(norm(A + U*V' - full(H)), '<', (norm(A) + norm(U*V')) * tol, ...
    'Sum of the hm representation of unstructured A and low-rank B');

H = hm('low-rank', U, V) + hm('banded', B, 3, 4);
CheckTestResult(norm(U*V' + B - full(H)), '<', (norm(A) + norm(U*V')) * tol, ...
    'Sum of the hm representation of low-rank A and banded B');

H = hm('low-rank', U, V) + hm('low-rank', W, Z);
CheckTestResult(norm(U*V' + W*Z' - full(H)), '<', (norm(U*V') + norm(W*Z')) * tol, ...
    'Sum of the hm representation of low-rank A and low-rank B');

A = randn(n, n); hmA = hm(A); 
B = randn(n, n); hmB = hm(B); hssB = hss(B);
H = hm(A) * hm(B);
CheckTestResult(norm(A * B - full(H)), '<', norm(A) * norm(B) * tol, ...
    'Product of the hm representation of unstructured A and unstructured B');

H = hmA * B;
CheckTestResult(norm(A * B - full(H)), '<', norm(A) * tol && isfloat(H), ...
    'Product of the HM representation of  A and dense B');

H = A * hmB;
CheckTestResult(norm(A * B - full(H)), '<', norm(A) * tol && isfloat(H), ...
    'Product of dense  A and HM representation of B');

H = hmA * hssB;
CheckTestResult(norm(A * B - full(H)), '<', norm(A) * tol && isa(H, 'hm'), ...
    'Product of the HM representation of  A and HSS representation of  B');

B = tril(triu(randn(n), -3),4);
H = hm(A) * hm('banded', B, 3, 4);
CheckTestResult(norm(A * B - full(H)), '<', norm(A) * norm(B) * tol, ...
    'Product of the hm representation of unstructured A and banded B');

A = tril(triu(randn(n), -6),3);
H = hm('banded', A, 6, 3) * hm('banded', B, 3, 4);
CheckTestResult(norm(A * B - full(H)), '<', norm(A) * norm(B) * tol, ...
    'Product of the hm representation of banded A and banded B');

U = randn(n, 6); V = randn(n, 6);
W = randn(n, 6); Z = randn(n, 6);
A = randn(n);
H = hm(A) * hm('low-rank', U, V);
CheckTestResult(norm(A * U*V' - full(H)), '<', norm(A) * norm(U) * norm(V) * tol, ...
    'Product of the hm representation of unstructured A and low-rank B');

H = hm('low-rank', U, V) * hm('banded', B, 3, 4);
CheckTestResult(norm(U*V' * B - full(H)), '<', norm(U) * norm(V) * norm(B) * tol, ...
    'Product of the hm representation of low-rank A and banded B');

H = hm('low-rank', U, V) * hm('low-rank', W, Z);
CheckTestResult(norm(U*V' * W*Z' - full(H)), '<', norm(U) * norm(V) * norm(W) * norm(Z) * tol, ...
    'Product of the hm representation of low-rank A and low-rank B');

% Matrix vector multiplication
H = hmgallery('rand', n, 10);
x = rand(n, 3); b = H * x;
CheckTestResult(norm(b - full(H)*x), '<', norm(H) * norm(x) * eps, ...
    'Matrix vector multiplication');

% Hadamard product of random hm matrices
A = hmgallery('rand', n, 10);
B = hmgallery('rand', n, 10);
H = A .* B;
CheckTestResult(norm(full(A) .* full(B) - full(H)), '<', norm(full(A) .* full(B)) * tol, ...
    'Hadamard product of the hm representation of random hm A and  B');

A = randn(n, n); hmA = hm(A); 
B = randn(n, n); hmB = hm(B); hssB = hss(B);
H = hmA .* B;
CheckTestResult(norm(A .* B - full(H)), '<', norm(A) * tol && isfloat(H), ...
    'Hadamard product of the HM representation of  A and dense B');

H = A .* hmB;
CheckTestResult(norm(A .* B - full(H)), '<', norm(A) * tol && isfloat(H), ...
    'Hadamard product of dense  A and HM representation of B');

H = hmA .* hssB;
CheckTestResult(norm(A .* B - full(H)), '<', norm(A) * tol && isa(H, 'hm'), ...
    'Hadamard product of the HM representation of  A and HSS representation of  B');

A = randn(n, n); hmA = hm(A);
CheckTestResult(norm(tril(A) - tril(hmA)), '<', norm(A) * tol, ...
    'Tril of the HM representation of A ');

CheckTestResult(norm(triu(A) - full(triu(hmA))), '<', norm(A) * tol, ...
    'Triu of the HM representation of A ');

% Power of an HM matrix
A = hmgallery('rand', n, 10);
H = A^2;
CheckTestResult(norm(full(A)^2 - full(H)), '<', norm(full(A))^2 * tol, ...
    'Square of the hm representation of random hm A');

A = hmgallery('rand', n, 10);
H = A^3;
CheckTestResult(norm(full(A)^3 - full(H)), '<', norm(full(A))^3 * tol, ...
    'Cube of the hm representation of random hm A');

for p = [ 4, 7, 9, 12, 15 ]
    A = hmgallery('rand', n, 10);
    A = A / norm(A); % Make sure things stay bounded
    H = A^p;
    CheckTestResult(norm(full(A)^p - full(H)), '<', norm(full(A))^p * tol, ...
        sprintf('%d-th Power of the hm representation of random hm A', p));
end

% Linear systems
H = hmgallery('rand', n, 3);
x = rand(n, 5); b = H*x; y = H \ b;
CheckTestResult(norm(x - y), '<', cond(full(H)) * norm(b), ...
    'Solution of a linear system (implicit LU)');

[L, U] = lu(H); y = U \ (L \ b);
CheckTestResult(norm(x - y), '<', cond(full(H)) * norm(b), ...
    'Solution of a linear system (explicit LU)');



end
