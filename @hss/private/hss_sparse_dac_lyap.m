function X = hss_sparse_dac_lyap(A, B, C, sA, sB)
% HSS_DAC_LYAP Divide and conquer method for solving A X + X B + C = 0 
%          where the coefficient matrices are represented both in the HSS format and in sparse format
kmax = 65;

debug = 0;
tol = 1e-8;

n = size(A, 1);

if A.leafnode == 1
	X = hss();
	X.D = lyap(A.D, B.D, C.D);

	X.topnode = 1;
	X.leafnode = 1;

	return;
end

X = blkdiag(...
	hss_sparse_dac_lyap(A.hssl, B.hssl, C.hssl, sA(1:A.ml, 1:A.nl), sB(1:A.ml, 1:A.nl)), ...
	hss_sparse_dac_lyap(A.hssr, B.hssr, C.hssr, sA(A.ml+1:end, A.nl+1:end), sB(A.ml+1:end, A.nl+1:end)) ...
);

[CU,CV] = hss_offdiag(C);
[AU,AV] = hss_offdiag(A);
[BU,BV] = hss_offdiag(B);

u = [ CU , AU , X * BU ];
v = [ CV , X' * AV, BV ];

[u, v] = compress_factors(u, v, 1.0);

A.topnode = 1;
B.topnode = 1;
%[LA,UA] = lu(sA);
%[LB,UB] = lu(sB);
%[Xu, Xv] = kpik_sylv(sA, LA, UA, sB, LB, UB, -u, v, 100, tol);
%[Xu, Xv] = kpik_sylv(A, speye(size(A)), A, B, speye(size(B)), B, -u, v, 100, tol);
[Xu, Xv] = SylvKrylov2(sA, sB, u, v, inf, tol);
% XX = lyap(full(A),full(B), -u*v');
% [Xu,D,Xv] = tsvd(XX,1e-12); Xu=Xu*D;
 %norm(full(sA * Xu * Xv' + Xu * (Xv' * sB') - u * v')) / norm(u * v')

A.topnode = 0;
B.topnode = 0;

X = X + hss('low-rank', Xu, Xv);
