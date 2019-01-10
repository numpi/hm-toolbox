function C = hss_sum(A, B, compress)
C = hss_sum_ric(A, B);
if ~exist('compress', 'var')
    compress = true;
end

if compress
    tol = hssoption('threshold');
    C = hss_compress(C,tol);
end

end

function C = hss_sum_ric(A,B)
C = A;
if  A.topnode ~= B.topnode || A.leafnode ~= B.leafnode
    error('HSS_SUM:: the two hss matrices have not compatible partitioning')
else
    % C.topnode = A.topnode;
    % C.leafnode = A.leafnode;
    if C.leafnode == 0
        if A.ml ~= B.ml || A.nl ~= B.nl ||A.mr ~= B.mr ||A.nr ~= B.nr
            error('HSS_SUM:: the two hss matrices have not compatible partitioning')
        end
        C.ml = A.ml;
        C.nl = A.nl;
        C.mr = A.mr;
        C.nr = A.nr;
    end
end

if C.leafnode == 1
    C.D = A.D + B.D;
    C.U = [A.U, B.U];
    C.V = [A.V, B.V];
else
    C.B12 = blkdiag(A.B12, B.B12);
    C.B21 = blkdiag(A.B21, B.B21);
    if C.topnode == 0
        C.Rl = blkdiag(A.Rl, B.Rl);
        C.Rr = blkdiag(A.Rr, B.Rr);
        C.Wl = blkdiag(A.Wl, B.Wl);
        C.Wr = blkdiag(A.Wr, B.Wr);
    end
    C.A11 = hss_sum_ric(A.A11, B.A11);
    C.A22 = hss_sum_ric(A.A22, B.A22);
end
end
