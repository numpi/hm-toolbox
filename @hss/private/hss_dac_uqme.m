function X = hss_dac_uqme(A, B, C, tol, debug, nrmtype)
% solve the UQME A X^2 + B X  + C = 0 with divide & conquer
% sparseA, sparseB additional sparse version of A and B ---> speed up in building krylov subspaces
if ~exist('tol', 'var')
	tol = hodlroption('threshold') * 1e2;
end
if ~exist('debug', 'var')
	debug = false;
end
if ~exist('nrmtype', 'var')
	nrmtype = 2;
end
[n, m] = size(A);
X = hss();

if ~isempty(A.D)	
	X.D = cr_uqme(A.D, B.D, C.D); 
	X.leafnode = 1;
	X.topnode = 1;
	if debug
		fprintf('Base case: Dimension: %d Residual = %e\n', n, norm(A.D * X.D^2 + B.D * X.D + C.D));
	end
else
	A.topnode = 1;
	B.topnode = 1;
	C.topnode = 1;
	X0 = blkdiag(hss_dac_uqme(A.A11, B.A11, C.A11, tol, debug, nrmtype),  hss_dac_uqme(A.A22, B.A22, C.A22, tol, debug, nrmtype));
	
	[AU, AV] = hss_offdiag(A, 'all');
	[BU, BV] = hss_offdiag(B, 'all');
	[CU, CV] = hss_offdiag(C, 'all');
	AV = X0' * (X0' * AV); 
	BV = X0' * BV;
	[u, v] = compress_factors([AU, BU, CU], [AV, BV, CV]);
	[dXu, dXv] = ek_uqme(A, B, X0, u, v, min(1000, n), tol, debug, nrmtype);
	if debug
		dX = hss('low-rank', dXu, dXv); 
		fprintf('Correction equation: Dimension %d, Rank of the sol. = %d, Residual = %e\n', n, size(dXu, 2), norm((A * dX  + A * X0 + B) * dX + A * dX * X0 + hss('low-rank', u, v)))
    end
	X = X0 + hss('low-rank', dXu, dXv);
	if debug 
		fprintf('Complete equation: Dimension %d, HSS rank of the sol.  = %d,  Res = %e\n', n, hssrank(X), norm((A * X  + B) * X + C))
	end	
end
