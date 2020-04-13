function [G,F] = cr_uqme(A2, A1, A0, tol, maxit)
% [G,F]=CR_UQME(A0,A1,A2) solves the UQMEs A0 + A1 X + A2 X^2 = 0 and 
% A2 + A1 Y + A0 Y^2 = 0 by means of cyclic reduction
% A0, A1, A2: matrix coefficients of the UQMEs
% G: solution of the UQME A0 + A1 X + A2 X^2 = 0
% F: solution of the UQME A2 + A1 Y + A0 Y^2 = 0
if ~exist('tol', 'var')
	tol = 1e-13;
end
if ~exist('maxit', 'var')
	maxit = 50;
end
n = size(A0, 1);
% initialization
AH = A1;
AT = A1;

B0 = A0; B1 = A1; B2 = A2;
% CR step
err = 1;
k = 0;
while err > tol && k < maxit
    F  = [B0; B2] / B1;
    F0 = F(1:n, :);
    F2 = F(n + 1:2 * n, :);
    W = F2 * B0;
    B0 = F0 * B0;
    AH = AH - W;
    B1 = B1 - W;
    W = F0 * B2;
    B2 = F2 * B2;
    AT = AT - W;

    B1 = B1 - W;
    err = min(norm(B0, 'fro'),norm(B2, 'fro'));
    k = k + 1;
end
G = -AH \ A0; 
F = -AT \ A2;

if k == maxit
    disp('CR_UQME:: Warning: reached the maximum number of iterations\n')
end
