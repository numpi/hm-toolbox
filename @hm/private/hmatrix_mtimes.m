function H = hmatrix_mtimes(H1, H2)
%HMATRIX_MTIMES Multiply two H matrices.

if size(H1, 2) ~= size(H2, 1)
    error('A * B: Dimension mismatch');
end

H = H1;

H.sz = [H1.sz(1), H2.sz(2)];

if ~isempty(H1.F)
    H.F = H1.F * H2.F;
else    
    H.A11 = hmatrix_mtimes(H1.A11, H2.A11);
    U = H1.U12 * (H1.V12.' * H2.U21);
    H.A11 = hmatrix_rank_update(H.A11, U, H2.V21);
    
    H.A22 = hmatrix_mtimes(H1.A22, H2.A22);
    U = H1.U21 * (H1.V21.' * H2.U12);
    H.A22 = hmatrix_rank_update(H.A22, U, H2.V12);
    
    H.U12 = [ hmatrix_mtimes_dense(H1.A11, H2.U12), H1.U12 ];
    H.V12 =	[ H2.V12, dense_mtimes_hmatrix(H1.V12.', H2.A22).' ];
    
    % [H.U12, H.V12] = compress_factors(H.U12, H.V12, norm(H.U12, 'fro') * norm(H.V12, 'fro'));
    
    H.U21 = [ hmatrix_mtimes_dense(H1.A22, H2.U21), H1.U21 ];
    H.V21 = [ H2.V21, dense_mtimes_hmatrix(H1.V21.', H2.A11).' ];
    
    % [H.U21, H.V21] = compress_factors(H.U21, H.V21, norm(H.U21, 'fro') * norm(H.V21, 'fro'));
end

end

