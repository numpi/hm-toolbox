function hss_TestOperations
clear

hssoption('block-size', 32);
hodlroption('block-size', 32);
n = 512;

tol = hssoption('threshold');
hnrm = @(A) norm(A, hssoption('norm'));
C = @(A) sqrt(max(size(A)));

A = randn(n);
H = hss(A)';
CheckTestResult(norm(A' - full(H)), '<', C(A) * hnrm(A) * tol, ...
    'Transposition for unstructured HSS');

A = tril(triu(randn(n), -4), 6);
H = hss('banded', A, 4, 6)';
CheckTestResult(norm(A' - full(H)), '<', C(A) * hnrm(A) * tol, ...
    'Transposition for banded HSS');

U = randn(n,7); V = randn(n,7);
H = hss('low-rank', U, V)';
CheckTestResult(norm(V * U' - full(H)), '<', C(A) * hnrm(U*V') * tol, ...
    'Transposition for low-rank HSS');

A = randn(n, n); hssA = hss(A);
B = randn(n, n); hssB = hss(B); hodlrB = hodlr(B);
H = hss(A) + hss(B);
CheckTestResult(norm(A + B - full(H)), '<', C(A) * (hnrm(A + B)) * tol, ...
    'Sum of the HSS representation of unstructured A and unstructured B');

H = hssA + B;
CheckTestResult(norm(A + B - full(H)), '<', C(A) * hnrm(A + B) * tol, ...
    'Sum of the HSS representation of  A and dense B');

H = A + hssB;
CheckTestResult(norm(A + B - full(H)), '<', C(A) * hnrm(A + B) * tol, ...
    'Sum of dense  A and HSS representation of B');

H = hssA + hodlrB;
CheckTestResult(norm(A + B - full(H)), '<', C(A) * hnrm(hssA + hssB) * tol, ...
    'Sum of the HSS representation of  A and hodlr representation of  B');

H = hssA - B;
CheckTestResult(norm(A - B - full(H)), '<', C(A) * hnrm(A - B) * tol, ...
    'Difference of the HSS representation of A and dense B');

H = A - hssB;
CheckTestResult(norm(A - B - full(H)), '<', C(A) * hnrm(A - B) * tol, ...
    'Difference of dense A and HSS representation of B');

H = hssA - hodlrB;
CheckTestResult(norm(A - B - full(H)), '<', C(A) * hnrm(A - B) * tol, ...
    'Difference of the HSS representation of A and hodlr representation of  B');

B = tril(triu(randn(n), -3),4);
H = hss(A) + hss('banded', B, 3, 4);
CheckTestResult(norm(A + B - full(H)), '<', C(A) * hnrm(B + A) * tol, ...
    'Sum of the HSS representation of unstructured A and banded B');

A = tril(triu(randn(n), -6),3);
H = hss('banded', A, 6, 3) + hss('banded', B, 3, 4);
CheckTestResult(norm(A + B - full(H)), '<', 2 * C(A) * hnrm(A + B) * tol, ...
    'Sum of the HSS representation of banded A and banded B');

U = randn(n, 6); V = randn(n, 6);
W = randn(n, 6); Z = randn(n, 6);
A = randn(n);
H = hss(A) + hss('low-rank', U, V);
CheckTestResult(norm(A + U*V' - full(H)), '<', C(A) * hnrm(A + U*V') * tol, ...
    'Sum of the HSS representation of unstructured A and low-rank B');

H = hss('low-rank', U, V) + hss('banded', B, 3, 4);
CheckTestResult(norm(U*V' + B - full(H)), '<', C(A) * hnrm(B + U*V') * tol, ...
    'Sum of the HSS representation of low-rank A and banded B');

H = hss('low-rank', U, V) + hss('low-rank', W, Z);
CheckTestResult(norm(U*V' + W*Z' - full(H)), '<', C(A) * hnrm(U*V' + W*Z') * tol, ...
    'Sum of the HSS representation of low-rank A and low-rank B');

A = randn(n, n); hssA = hss(A);
B = randn(n, n); hssB = hss(B); hodlrB = hodlr(B);
H = hss(A) * hss(B);
CheckTestResult(norm(A * B - full(H)), '<', C(A) * hnrm(A * B) * tol, ...
    'Product of the HSS representation of unstructured A and unstructured B');

H = hssA * B;
CheckTestResult(norm(A * B - full(H)), '<', C(A) * hnrm(A * B) * tol, ...
    'Product of the HSS representation of  A and dense B');

H = A * hssB;
CheckTestResult(norm(A * B - full(H)), '<', C(A) * hnrm(A * B) * tol, ...
    'Product of dense  A and HSS representation of B');

H = hssA * hodlrB;
CheckTestResult(norm(A * B - full(H)), '<', C(A) * hnrm(A * B) * tol, ...
    'Product of the HSS representation of  A and hodlr representation of  B');

B = tril(triu(randn(n), -3),4);
H = hss(A) * hss('banded', B, 3, 4);
CheckTestResult(norm(A * B - full(H)), '<', C(A) * hnrm(A * B) * tol, ...
    'Product of the HSS representation of unstructured A and banded B');

A = tril(triu(randn(n), -6),3);
H = hss('banded', A, 6, 3) * hss('banded', B, 3, 4);
CheckTestResult(norm(A * B - full(H)), '<', C(A) * hnrm(A * B) * tol, ...
    'Product of the HSS representation of banded A and banded B');

U = randn(n, 6); V = randn(n, 6);
W = randn(n, 6); Z = randn(n, 6);
A = randn(n);
H = hss(A) * hss('low-rank', U, V);
CheckTestResult(norm(A * U*V' - full(H)), '<', hnrm(A * U * V') * tol, ...
    'Product of the HSS representation of unstructured A and low-rank B');

H = hss('low-rank', U, V) * hss('banded', B, 3, 4);
CheckTestResult(norm(U*V' * B - full(H)), '<', hnrm(U * V' * B) * tol, ...
    'Product of the HSS representation of low-rank A and banded B');

H = hss('low-rank', U, V) * hss('low-rank', W, Z);
CheckTestResult(norm(U*V' * W*Z' - full(H)), '<', hnrm(U * V'* W * Z') * tol, ...
    'Product of the HSS representation of low-rank A and low-rank B');

% Matrix vector multiplication
H = hssgallery('rand', n, 10);
x = rand(n, 3); b = H * x;
CheckTestResult(norm(b - full(H)*x), '<', hnrm(H) * norm(x) * eps, ...
    'Matrix vector multiplication');

% Hadamard product of random hss matrices
A = hssgallery('rand', n, 10);
B = hssgallery('rand', n, 10);
H = A .* B;
CheckTestResult(norm(full(A) .* full(B) - full(H)), '<', hnrm(full(A) .* full(B)) * tol, ...
    'Hadamard product of the HSS representation of random hss A and  B');

A = randn(n, n); hssA = hss(A);
B = randn(n, n); hssB = hss(B); hodlrB = hodlr(B);
H = hssA .* B;
CheckTestResult(norm(A .* B - full(H)), '<', hnrm(A .* B) * tol, ...
    'Hadamard product of the HSS representation of  A and dense B');

H = A .* hssB;
CheckTestResult(norm(A .* B - full(H)), '<', hnrm(A .* B) * tol, ...
    'Hadamard product of dense  A and HSS representation of B');

H = hssA .* hodlrB;
CheckTestResult(norm(A .* B - full(H)), '<', hnrm(A .* B) * tol, ...
    'Hadamard product of the HSS representation of  A and hodlr representation of  B');

A = randn(n, n); hssA = hss(A);
CheckTestResult(norm(tril(A) - tril(hssA)), '<', hnrm(A) * tol, ...
    'Tril of the HSS representation of A ');

CheckTestResult(norm(triu(A) - full(triu(hssA))), '<', hnrm(A) * tol, ...
    'Triu of the HSS representation of A ');

% Power of an hodlr matrix
A = hssgallery('rand', n, 10);
H = A^2;
CheckTestResult(norm(full(A)^2 - full(H)), '<', hnrm(full(A))^2 * tol, ...
    'Square of the hss representation of random hss A');

A = hssgallery('rand', n, 10);
H = A^3;
CheckTestResult(norm(full(A)^3 - full(H)), '<', hnrm(full(A))^3 * tol, ...
    'Cube of the hss representation of random hss A');

for p = [ 4, 7, 9, 12, 15 ]
    A = hssgallery('rand', n, 10);
    A = A / norm(A); % Make sure things stay bounded
    H = A^p;
    CheckTestResult(norm(full(A)^p - full(H)), '<', hnrm(full(A))^p * tol, ...
        sprintf('%d-th Power of the hss representation of random hss A', p));
end

% Linear systems
H = hssgallery('rand', n, 3);
x = rand(n, 5); b = H*x; y = H \ b;
CheckTestResult(norm(x - y), '<', cond(full(H)) * norm(b) * tol, ...
    'Solution of a linear system (implicit ULV)');

F = ulv(H); y = ulv_solve(F, b);
CheckTestResult(norm(x - y), '<', cond(full(H)) * norm(b) * tol, ...
    'Solution of a linear system (explicit ULV)');

% Inversion and backslash
B = hssgallery('rand', n, 3);
X = H \ B;
CheckTestResult(norm(H*X - B), '<', (hnrm(X)*hnrm(H) + hnrm(B)) * tol, ...
    'Backslash operator for HSS A \ B (hss rank 3)');

% Inversion and backslash
B = hssgallery('rand', n, 3);
X = B / H;
CheckTestResult(norm(X*H - B), '<', (hnrm(X)*hnrm(H) + hnrm(B)) * tol, ...
    'Backslash operator for HSS A / B (hss rank 3)');

X = inv(H);
CheckTestResult(norm(H*X - eye(n, 'like', H)), '<', (hnrm(X)*hnrm(H) + hnrm(B)) * tol, ...
    'Inversion (hss rank 3)');

H = hssgallery('rand', n, 12);
B = hssgallery('rand', n, 12);
X = H \ B;
CheckTestResult(norm(H*X - B), '<', (hnrm(X)*hnrm(H) + hnrm(B)) * tol, ...
    'Backslash operator for HSS A \ B (hss rank 12)');

X = inv(H);
CheckTestResult(norm(H*X - eye(n, 'like', H)), '<', (hnrm(X)*hnrm(H) + hnrm(B)) * tol, ...
    'Inversion (hss rank 12)');

% Inversion with non-standard clustering
c = [ 21  42  42  42  42  42  42  42  59  76  93 111 128 145 162 180 ];
n = c(end);
H = hss(rand(n), 'cluster', c);
X = inv(H);
CheckTestResult(norm(H*X - eye(n, 'like', H)), '<', (hnrm(X)*hnrm(H) + hnrm(B)) * tol, ...
    'Inversion (non-standard cluster)');

% Random cluster
c = cumsum(randi(32, 1, 16));
n = c(end);
H = hss(rand(n), 'cluster', c);
X = inv(H);
CheckTestResult(norm(H*X - eye(n, 'like', H)), '<', (hnrm(X)*hnrm(H) + hnrm(B)) * tol, ...
    'Inversion (non-standard random cluster)');


end
