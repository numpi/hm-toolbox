function H = hodlr_mtimes(H1, H2, compress)
%HODLRATRIX_MTIMES Multiply two H matrices.
if ~exist('compress', 'var')
    compress = true;
end
if size(H1, 2) ~= size(H2, 1)
    error('A * B: Dimension mismatch');
end
nrm = normest_Afun(@(x) H1 * (H2 * x), @(x) H2' * (H1' * x), size(H2, 2)); % estimate the norm of the product
H = hodlr_mtimes_ric(H1, H2, nrm);

if compress
    H = compress_hodlr(H, nrm);
end
end

function H = hodlr_mtimes_ric(H1, H2, nrm)
H = H1;

H.sz = [H1.sz(1), H2.sz(2)];

if is_leafnode(H1)
    H.F = H1.F * H2.F;
else    
    H.A11 = hodlr_mtimes_ric(H1.A11, H2.A11, nrm);
    U = H1.U12 * (H1.V12.' * H2.U21);
    H.A11 = hodlr_rank_update(H.A11, U, H2.V21, nrm);
    
    H.A22 = hodlr_mtimes_ric(H1.A22, H2.A22, nrm);
    U = H1.U21 * (H1.V21.' * H2.U12);
    H.A22 = hodlr_rank_update(H.A22, U, H2.V12, nrm);
    
    H.U12 = [ hodlr_mtimes_dense(H1.A11, H2.U12), H1.U12 ];
    H.V12 = [ H2.V12, dense_mtimes_hodlr(H1.V12.', H2.A22).' ];
    
    % [H.U12, H.V12] = compress_factors(H.U12, H.V12, norm(H.U12, 'fro') * norm(H.V12, 'fro'));
    
    H.U21 = [ hodlr_mtimes_dense(H1.A22, H2.U21), H1.U21 ];
    H.V21 = [ H2.V21, dense_mtimes_hodlr(H1.V21.', H2.A11).' ];
    
    % [H.U21, H.V21] = compress_factors(H.U21, H.V21, norm(H.U21, 'fro') * norm(H.V21, 'fro'));
end

end

