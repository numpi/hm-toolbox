function [X, converged] = newton_nare(A, B, C, D, X0, debug)
% X=NEWTON_NARE(A,B,C,D,X0) solves the NARE C + XA + DX - XBX = 0
% by means of Newton's method 
%    A, B, C, D: matrix coefficients
%    X0: initial approximation
%    X : solution of the NARE
if ~exist('debug', 'var')
	debug = false;
end
tol = 1e-13;
kmax = 20;
X = X0;
err = 1; err_old = err;
k = 1;
converged = true;
linesearch = 1;
RX = C + X * A + D * X - X * B * X;
while err > tol && k < kmax 
    H = lyap(D - X * B, A - B * X, RX);
    if linesearch
    	V = H*B*H;
    	a = trace(RX' * RX);
    	b = trace(RX' * V);
    	c = trace(V' * V);
    	tk = fminbnd(@(t) a*(1-t)^2-2*b*(1-t)*t^2+c*t^4,0,2);
    else
	tk = 1;
    end
    X = X + tk * H;
    RX = C + X * A + D * X - X * B * X;
    err_old = min(err, err_old);
    err = norm(RX, 'fro')/norm(X, 'fro');
    if  err_old < err && k > 2 
	if err > 1e-10
 		converged = false;
		%fprintf('NEWTON_NARE::Warning: non descending sequence, Dim = %d, It = %d,  Rel. Res = %e, Abs Res = %e\n', size(A, 1), k, err, err * norm(X,'fro'))
	end
	return;
    end
    k = k + 1;
end
if k == kmax 
	if err > 1e-10
    		converged = false;
	end
	if debug
    		fprintf('NEWTON_NARE::Warning: reached the maximum number of iterations Res = %e\n', err)
	end
end
