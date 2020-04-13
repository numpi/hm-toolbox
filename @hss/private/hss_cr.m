function [G,F] = hss_cr(A2, A1, A0, tol, debug, nrmtype, maxit)
% [G,F] = HSS_CR(A0,A1,A2) solves the UQMEs A0 + A1 X + A2 X^2 = 0 and 
% A2 + A1 Y + A0 Y^2 = 0 by means of cyclic reduction
% A0, A1, A2: matrix coefficients of the UQMEs
% G: solution of the UQME A0 + A1 X + A2 X^2 = 0
% F: solution of the UQME A2 + A1 Y + A0 Y^2 = 0
if ~exist('tol', 'var')
	tol = 1e-13;
end
if ~exist('debug', 'var')
	debug = 0;
end
if ~exist('maxit', 'var')
	maxit = 30;
end
if ~exist('nrmtype', 'var')
	nrmtype = 2;
end
n = size(A0,1);
% initialization
AH = A1;
if nargout == 2
	AT = A1;
end
B0 = A0; B1 = A1; B2 = A2;
% CR step
err = 1;
k = 0;

while err > tol && k < maxit
    F0 = B0 / B1;
    F2 = B2 / B1;
    W = F2 * B0;
    B0 = F0 * B0;
    AH = AH - W;
    B1 = B1 - W;
    W = F0 * B2;
    B2 = F2 * B2;
if nargout == 2
    AT = AT - W;
end
    B1 = B1 - W;
    %err = min(norm(B0),norm(B2));
    %G = -AH \ A0;
    %err = (A2 * G + A1) * G + A0;
    %err = norm(err, nrmtype); 
    err = min(norm(B0, 'fro'), norm(B2, 'fro'));
    if debug
	fprintf('CR iteration: %d, residue: %e\n', k, err)
    end
    k = k + 1;
end
G = -AH \ A0;
if nargout == 2
	F = -AT \ A2;
end

if k == maxit
    disp('HSS_CR:: Warning: reached the maximum number of iterations')
end
