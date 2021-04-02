function H = halr_rank_update(H, U, V, nrm)
%HALR_RANK_UPDATE Perform a low rank update to H.
%
% H = HALR_RANK_UPDATE(H, U, V) computes an HODLR representation of the
%     matrix H + U * V'; Truncation is performed on each off-diagonal block
%     with respect to its 2-norm.
%
% H = HALR_RANK_UPDATE(H, U, V, NRM) truncates the result of H + U*V'
%     with respect to the global threshold NRM * HODLROPTION('threshold');

if is_leafnode(H)
    if H.admissible
        if exist('nrm', 'var')
            if isempty(nrm)
                [H.U, H.V] = compress_factors([H.U, U], [H.V, V]);
            else
                [H.U, H.V] = compress_factors([H.U, U], [H.V, V], nrm);
            end
        else
            [H.U, H.V] = compress_factors([H.U, U], [H.V, V]);
        end
    else
        H.F = H.F + U * V';
    end
else
    m1 = H.A11.sz(1);
    n1 = H.A11.sz(2);
    if exist('nrm', 'var')
        H.A11 = halr_rank_update(H.A11, U(1:m1,:), V(1:n1,:), nrm);
        H.A12 = halr_rank_update(H.A12, U(1:m1,:), V(n1+1:end,:), nrm);
        H.A21 = halr_rank_update(H.A21, U(m1+1:end,:), V(1:n1,:), nrm);
        H.A22 = halr_rank_update(H.A22, U(m1+1:end,:), V(n1+1:end,:), nrm);
    else
        H.A11 = halr_rank_update(H.A11, U(1:m1,:), V(1:n1,:));
        H.A12 = halr_rank_update(H.A12, U(1:m1,:), V(n1+1:end,:));
        H.A21 = halr_rank_update(H.A21, U(m1+1:end,:), V(1:n1,:));
        H.A22 = halr_rank_update(H.A22, U(m1+1:end,:), V(n1+1:end,:));
    end
end

