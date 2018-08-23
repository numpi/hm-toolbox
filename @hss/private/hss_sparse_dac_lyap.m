function X = hss_sparse_dac_lyap(A, C, sA, use_sylv)
% HSS_DAC_LYAP Divide and conquer method for solving A X + X B + C = 0
%          where the coefficient matrices are represented both in the HSS format and in sparse format

diagblk = false; % enables the trick suggested by Beckermann for the correction equation

if A.leafnode == 1
    X = hss();
    X.D = lyap(A.D, C.D);
    
    X.topnode = 1;
    X.leafnode = 1;
    
    return;
end

if ~exist('use_sylv', 'var')
    use_sylv = true;
end

X = blkdiag(...
    hss_sparse_dac_lyap(A.A11, C.A11, sA(1:A.ml, 1:A.nl), use_sylv), ...
    hss_sparse_dac_lyap(A.A22, C.A22, sA(A.ml+1:end, A.nl+1:end), use_sylv) ...
    );

[CU, CV] = hss_offdiag(C);
[AU, AV] = hss_offdiag(A);

u = [ CU , AU , X * AV ];
v = [ CV , X' * AV, AU ];

tol = hssoption('threshold');

if use_sylv
    if issymmetric(sA)
        AA = ek_struct(sA, true);
        AT = AA;
    else
        [AA, AT] = ek_struct(sA);
    end
    
    [u, v] = compress_factors(u, v, 1.0);
    
    if size(u,2) > 0
        if ~diagblk
            [Xu,Xv] = ek_sylv(AA, AA, u, v, inf, tol);
        else
            [Xu,Xv] = ek_sylv_blkdiag2(sA, sA, u, v, inf, tol);
        end
        A.topnode = 0;
        X = X + hss('low-rank', Xu, Xv);
    end
else
    % Find a symmetric definite decomposition u*v'= Up*Up' - Un*Un'
    [Up, Un] = ek_definite_splitting(u, v, tol);
    
    A.topnode = 1;
    
    if ~isempty(Up)
        [Xup, ~, AA] = ek_lyap(sA, Up, inf, tol);
    else
        Xup = Up;
    end
    
    if ~isempty(Un)
        if exist('AA', 'var')
            Xun = ek_lyap(AA, Un, inf, tol);
        else
            Xun = ek_lyap(sA, Un, inf, tol);
        end
    else
        Xun = Un;
    end
    
    Xu = [ Xup,  Xun ];
    Xv = [ Xup, -Xun ];
    A.topnode = 0;
    X = X + hss('low-rank', Xu, Xv);
end



