function H = hmatrix_repartition(H)
%REPARTITION_HMATRIX 

if is_leafnode(H)
    if ~H.admissible
        H = hmatrix('adaptive', H.F, size(H.F, 1), size(H.F, 2));
        %[U, S, V] = tsvd(H.F, hmatrixoption('threshold'));
        %if size(U, 2) < min(size(H.F) / 2)
        %    H.admissible = true;
        %    H.U = U * S;
        %    H.V = V;
        %end
    end
else
    H.A11 = hmatrix_repartition(H.A11);
    H.A12 = hmatrix_repartition(H.A12);
    H.A21 = hmatrix_repartition(H.A21);
    H.A22 = hmatrix_repartition(H.A22);
    
    if H.A11.admissible && H.A12.admissible && H.A21.admissible && H.A22.admissible
        [U, V] = merge_low_rank(H.A11, H.A12, H.A21, H.A22);
        
        mem = sum(H.sz) * size(U, 2);
        mem11 = sum(H.A11.sz) * size(H.A11.U, 2);
        mem12 = sum(H.A12.sz) * size(H.A12.U, 2);
        mem21 = sum(H.A21.sz) * size(H.A21.U, 2);
        mem22 = sum(H.A22.sz) * size(H.A22.U, 2);
        
        if mem <= mem11 + mem12 + mem21 + mem22
            H.A11 = []; H.A12 = []; H.A21 = []; H.A22 = [];
            H.U = U;
            H.V = V;
            H.admissible = true;
        end
    end
end

end

function [U, V] = merge_low_rank(A11, A12, A21, A22)
    [U, V] = compress_factors(blkdiag([ A11.U, A12.U ], [ A21.U, A22.U ]), ...
                [ blkdiag(A11.V, A12.V), blkdiag(A21.V, A22.V) ]);
end

