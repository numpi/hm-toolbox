function X = sparse_dac_lyap(A, B, C, sA, sB)
% LYAP_DAQ Divide and conquer method for solving A X + X B + C = 0
%          where all the matrices are represented in the HODLR format
%k = 3;
X = A;
n = size(A,1);
debug = 0;

if is_leafnode(A) && is_leafnode(B) && is_leafnode(C)
    X.F = lyap(A.F, B.F, C.F);
    return
end

if n <= 256
    X = hm(lyap(full(A), full(B), full(C)));
    return
end

X.A11 = sparse_dac_lyap(A.A11, B.A11, C.A11, sA(1:n/2,1:n/2),sB(1:n/2,1:n/2)); % Recursive solution on the diagonal blocks
X.A22 = sparse_dac_lyap(A.A22, B.A22, C.A22, sA(n/2+1:end,n/2+1:end),sB(n/2+1:end,n/2+1:end));
X.U21 = zeros(size(A.U21,1),0);
X.U12 = zeros(size(A.U12,1),0);
X.V21 = zeros(size(A.V21,1),0);
X.V12 = zeros(size(A.V12,1),0);
if debug
    AA=A; AA.U21 = zeros(size(A.U21,1),0); AA.U12 = zeros(size(A.U12,1),0); AA.V21 = zeros(size(A.V21,1),0); AA.V12 = zeros(size(A.V12,1),0);
    BB=B; BB.U21 = zeros(size(B.U21,1),0); BB.U12 = zeros(size(B.U12,1),0); BB.V21 = zeros(size(B.V21,1),0); BB.V12 = zeros(size(B.V12,1),0);
    CC=C; CC.U21 = zeros(size(C.U21,1),0); CC.U12 = zeros(size(C.U12,1),0); CC.V21 = zeros(size(C.V21,1),0); CC.V12 = zeros(size(C.V12,1),0);
    norm(AA*X+X*BB+CC,'fro')/norm(X,'fro')/norm(A,'fro')/sqrt(size(AA,1))
end

u1 = [C.U21, A.U21, X.A22 * B.U21];	% Compute low rank factorization of the right hand-side
v1 = [C.V21, X.A11' * A.V21, B.V21];
[u1, v1] = compress_factors(u1, v1, 1);%lr_norm(u1, v1)

u2 = [C.U12, A.U12, X.A11 * B.U12];
v2 = [C.V12, X.A22' * A.V12, B.V12];
[u2, v2] = compress_factors(u2, v2, 1);%lr_norm(u2, v2)
if size(u1,2) > 0 && size(u2,2) > 0
    u = [zeros(size(u2,1), size(u1, 2)), u2; u1, zeros(size(u1,1), size(u2, 2))];
    v = [ v1, zeros(size(v1,1), size(v2, 2));zeros(size(v2,1), size(v1, 2)), v2];
    if debug
        dA = A; dA.A11=hm(zeros(size(A.A11))); dA.A22=hm(zeros(size(A.A22))); dB = B; dB.A11=hm(zeros(size(B.A11))); dB.A22=hm(zeros(size(B.A22)));
        dC = C; dC.A11=hm(zeros(size(C.A11))); dC.A22=hm(zeros(size(C.A22)));
        norm(dC + dA * X + X * dB + hm('low-rank',u,v))
    end
    
    % Solve with Krylov methods for the low-rank update
    tol = hmoption('threshold');
    [ Xu, Xv ] = ek_sylv(sA, sB, u, v, inf, tol);
    if debug
        XX=hm('low-rank',Xu,Xv);
        norm(A*XX+XX*B+hm('low-rank',u,v),'fro')/norm(XX,'fro')/norm(A,'fro')/sqrt(size(XX,1))
    end
    X = hmatrix_rank_update(X, Xu, Xv);
end



