function H = hmatrix_rank_update(H, U, V, nrm)
%HMATRIX_RANK_UPDATE Perform a low rank update to H.
%
% H = HMATRIX_RANK_UPDATE(H, U, V) computes an HM representation of the
%     matrix H + U * V'; Truncation is performed on each off-diagonal block
%     with respect to its 2-norm.
%
% H = HMATRIX_RANK_UPDATE(H, U, V, NRM) truncates the result of H + U*V'
%     with respect to the global threshold NRM * HMOPTION('threshold'); 

if ~isempty(H.F)
    H.F = H.F + U * V.';
else
    m1 = H.A11.sz(1);
    n1 = H.A11.sz(2);
    
    if exist('nrm', 'var')
        [H.U12, H.V12] = compress_factors([ H.U12, U(1:m1,:) ], ...
            [ H.V12, V(n1+1:end,:) ], nrm);
        [H.U21, H.V21] = compress_factors([ H.U21, U(m1+1:end,:) ], ...
            [ H.V21, V(1:n1,:) ], nrm);
    else
        [H.U12, H.V12] = compress_factors([ H.U12, U(1:m1,:) ], ...
            [ H.V12, V(n1+1:end,:) ]);
        [H.U21, H.V21] = compress_factors([ H.U21, U(m1+1:end,:) ], ...
            [ H.V21, V(1:n1,:) ]);
    end
    
    H.A11 = hmatrix_rank_update(H.A11, U(1:m1,:), V(1:n1,:));
    H.A22 = hmatrix_rank_update(H.A22, U(m1+1:end,:), V(n1+1:end,:));
end

