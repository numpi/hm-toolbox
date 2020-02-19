function X = dac_lyap_blr(varargin)

A = varargin{1};

is_sparse = false;
is_lyapunov   = false;

switch nargin
    case 2
        C = varargin{2};
        is_lyapunov = true;
    case 3
        if issparse(varargin{3})
            sA = varargin{3};
            C = varargin{2};
            is_sparse = true;
            is_lyapunov = true;
        else
            B = varargin{2};
            C = varargin{3};
        end
    case 5
        B = varargin{2};
        C = varargin{3};
        sA = varargin{4};
        sB = varargin{5};
        is_sparse = true;
end

X = hmatrix;
X.sz = C.sz;

tol = hmatrixoption('threshold');

if is_leafnode(C) && C.admissible
    if is_sparse
        if is_lyapunov
            M = ek_struct(sA, issymmetric(sA));
            [X.U, X.V] = ek_sylv(M, M, C.U, C.V, inf, tol);
        else
            [X.U, X.V] = ek_sylv(sA, sB, C.U, C.V, inf, tol);
        end
    else
        if is_lyapunov
            M = ek_struct(A);
            [X.U, X.V] = ek_sylv(M, M, C.U, C.V, inf, tol);
        else
            [X.U, X.V] = ek_sylv(A, B, C.U, C.V, inf, tol);
        end
    end
    
    X.admissible = true;
elseif is_leafnode(A) && (is_lyapunov || is_leafnode(B))
    if is_lyapunov
        X.F = lyap(A.F, full(C));
    else
        X.F = lyap(A.F, B.F, full(C));
    end
    X.admissible = false;
elseif ~is_leafnode(A) && (is_lyapunov || ~is_leafnode(B))
            
    m1 = A.A11.sz(1); n1 = A.A11.sz(2);
    m2 = size(A) - m1; n2 = size(A, 2) - n1;
    
    if is_leafnode(C)       
        % Create an @hmatrix representation of C according to the splitting
        % of A and B, so we can perform the low-rank update approach
        % recursively. 
        CC = hmatrix; CC.sz = [ size(C, 1) , size(C, 2) ];
        CC.admissible = false;
                
        CC.A11 = hmatrix; CC.A11.sz = [ m1 n1 ]; CC.A11.F = C.F(1:m1, 1:n1);
        CC.A12 = hmatrix; CC.A12.sz = [ m1 n2 ]; CC.A12.F = C.F(1:m1, n1+1:end);
        CC.A21 = hmatrix; CC.A21.sz = [ m2 n1 ]; CC.A21.F = C.F(m1+1:end, 1:n1);
        CC.A22 = hmatrix; CC.A22.sz = [ m2 n2 ]; CC.A22.F = C.F(m1+1:end, n1+1:end);        
        
        C = CC;
    end
    
    if is_lyapunov && ~is_sparse
        X.A11 = dac_lyap_blr(A.A11, C.A11);
        X.A12 = dac_lyap_blr(A.A11, A.A22', C.A12);
        X.A21 = dac_lyap_blr(A.A22, A.A11', C.A21);
        X.A22 = dac_lyap_blr(A.A22, C.A22);
    elseif is_lyapunov && is_sparse
        X.A11 = dac_lyap_blr(A.A11, C.A11, sA(1:m1, 1:n1));
        X.A12 = dac_lyap_blr(A.A11, A.A22', C.A12, sA(1:m1, 1:n1), sA(m1+1:end, n1+1:end)');
        X.A21 = dac_lyap_blr(A.A22, A.A11', C.A21, sA(m1+1:end, n1+1:end), sA(1:m1, 1:n1)');
        X.A22 = dac_lyap_blr(A.A22, C.A22, sA(m1+1:end, n1+1:end));
    elseif ~is_lyapunov && ~is_sparse
        X.A11 = dac_lyap_blr(A.A11, B.A11, C.A11);
        X.A12 = dac_lyap_blr(A.A11, B.A22, C.A12);
        X.A21 = dac_lyap_blr(A.A22, B.A11, C.A21);
        X.A22 = dac_lyap_blr(A.A22, B.A22, C.A22);
    else
        X.A11 = dac_lyap_blr(A.A11, B.A11, C.A11, sA(1:m1, 1:n1), sB(1:m1, 1:n1));
        X.A12 = dac_lyap_blr(A.A11, B.A22, C.A12, sA(1:m1, 1:n1), sB(m1+1:end, n1+1:end));
        X.A21 = dac_lyap_blr(A.A22, B.A11, C.A21, sA(m1+1:end, n1+1:end), sB(1:m1, 1:n1));
        X.A22 = dac_lyap_blr(A.A22, B.A22, C.A22, sA(m1+1:end, n1+1:end), sB(m1+1:end, n1+1:end));
    end
    
    UA = -blkdiag(A.U12, A.U21); 
    VA = [ X.A21' * A.V12 , X.A11' * A.V21 ; ...
           X.A22' * A.V12 , X.A12' * A.V21 ];
       
    if is_lyapunov
        UB = [ X.A12 * A.V12 , X.A11 * A.V21 ; ...
           X.A22 * A.V12 , X.A21 * A.V21 ];
        VB = -blkdiag(A.U12, A.U21);
    else       
        UB = [ X.A12 * B.U21 , X.A11 * B.U12 ; ...
               X.A22 * B.U21 , X.A21 * B.U12 ];
        VB = -blkdiag(B.V21, B.V12);
    end
    
    [U, V] = compress_factors([ UA , UB ], [ VA, VB ]);
    
    if is_sparse
        if is_lyapunov
            M = ek_struct(sA);
            [Xu, Xv] = ek_sylv(M, M, -U, V, inf, tol);        
        else                    
            [Xu, Xv] = ek_sylv(sA, sB, -U, V, inf, tol);
        end
    else
        if is_lyapunov
            M = ek_struct(A);
            [Xu, Xv] = ek_sylv(M, M, -U, V, inf, tol);
        else                    
            [Xu, Xv] = ek_sylv(A, B, -U, V, inf, tol);
        end
    end
    
    X = rank_update(X, Xu, Xv);
else
    error('Non-conformal partitioning in A and B');        
end

end