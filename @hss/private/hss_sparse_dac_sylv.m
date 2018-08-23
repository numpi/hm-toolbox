function X = hss_sparse_dac_sylv(A, B, C, sA, sB)
% HSS_DAC_LYAP Divide and conquer method for solving A X + X B + C = 0
%          where the coefficient matrices are represented both in the HSS format and in sparse format

if A.leafnode == 1
    X = hss();
    
    X.D = lyap(A.D, B.D, C.D);
    
    X.topnode = 1;
    X.leafnode = 1;
    
    return;
end

X = blkdiag(...
    hss_sparse_dac_sylv(A.A11, B.A11, C.A11, sA(1:A.ml, 1:A.nl), sB(1:A.ml, 1:A.nl)), ...
    hss_sparse_dac_sylv(A.A22, B.A22, C.A22, sA(A.ml+1:end, A.nl+1:end), sB(A.ml+1:end, A.nl+1:end)) ...
    );

[CU,CV] = hss_offdiag(C);
[AU,AV] = hss_offdiag(A);
[BU,BV] = hss_offdiag(B);

u = [ CU , AU , X * BU ];
v = [ CV , X' * AV, BV ];

[u, v] = compress_factors(u, v, 1.0);

A.topnode = 1;
B.topnode = 1;

% [~,ru] = qr(u,0); [~,rv] = qr(v,0);
% tol = hssoption('threshold') / norm(ru * rv');
tol = hssoption('threshold');
[Xu, Xv] = ek_sylv(sA, sB, u, v, inf, tol);

A.topnode = 0;
B.topnode = 0;

X = X + hss('low-rank', Xu, Xv);
