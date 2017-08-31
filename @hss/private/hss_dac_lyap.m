function X = hss_dac_lyap(A,B,C)
% LYAP_DAQ Divide and conquer method for solving A X + X B + C = 0 
%          where all the matrices are represented in the HODLR format
k = 2;

debug = 0;
tol = 1e-8;

if A.leafnode == 1
	X = hss();
	X.D = lyap(A.D, B.D, C.D);

	X.topnode = 1;
	X.leafnode = 1;

	return;
end

X = blkdiag(...
	hss_dac_lyap(A.hssl, B.hssl, C.hssl), ...
	hss_dac_lyap(A.hssr, B.hssr, C.hssr) ...
);

[CU,CV] = hss_offdiag(C);
[AU,AV] = hss_offdiag(A);
[BU,BV] = hss_offdiag(B);

u = [ CU , AU , X * BU ];
v = [ CV , X' * AV, BV ];

[u, v] = compress_factors(u, v, norm(u) * norm(v));

A.topnode = 1;
B.topnode = 1;
[LA,UA] = lu(full(A));
[LLB,UB] = lu(full(B));
% [  Xu, Xv ] = kpik_sylv(A, speye(size(A)), A, B, speye(size(B)), B, u, v, 100, tol, tol);
[  Xu, Xv ] = kpik_sylv(full(A), LA, UA, full(B), LLB, UB, u, v, 100, tol, tol);
[ Xu, Xv] = kpik_sylv(-A, -LA, UA, -B, -LLB, UB, u, v, 100, tol, tol);
% [Xu, Xv] = SylvKrylov(A, B, u, v, 20);
norm(A * Xu * Xv' + Xu * (Xv' * B') - u * v') / norm(u * v')

A.topnode = 0;
B.topnode = 0;

X = X - hss('low-rank', Xu, Xv);
