function X = hss_dac_lyap(A,B,C)
% HSS_DAC_LYAP Divide and conquer method for solving A X + X B + C = 0
%          where all the matrices are represented in the HSS format

tol = hssoption('threshold');

if A.leafnode == 1
    X = hss();
    X.D = lyap(A.D, B.D, C.D);
    
    X.topnode = 1;
    X.leafnode = 1;
    
    return;
end

X = blkdiag(...
    hss_dac_lyap(A.A11, B.A11, C.A11), ...
    hss_dac_lyap(A.A22, B.A22, C.A22) ...
    );


[CU,CV] = hss_offdiag(C);
[AU,AV] = hss_offdiag(A);
[BU,BV] = hss_offdiag(B);


u = [ CU , AU , X * BU ];
v = [ CV , X' * AV, BV ];

[u, v] = compress_factors(u, v, norm(u) * norm(v));

A.topnode = 1;
B.topnode = 1;

[Xu, Xv] = ek_sylv(A, B, u, v, inf, tol);

A.topnode = 0;
B.topnode = 0;

X = X + hss('low-rank', Xu, Xv);
