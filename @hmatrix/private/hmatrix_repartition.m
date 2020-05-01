function H = hmatrix_repartition(H, maxrank, nrm)
%REPARTITION_HMATRIX

if ~exist('nrm', 'var')
    nrm = [];
end

if is_leafnode(H)
    if ~H.admissible
        if ~isempty(nrm)
            H = hmatrix_from_adaptive(H.F, size(H.F,1), size(H.F,2), maxrank, @(i,j) 1, nrm);
        else
            H = hmatrix_from_adaptive(H.F, size(H.F,1), size(H.F,2), maxrank);
        end
    else
        if size(H.U, 2) > min(maxrank, min(size(H.U, 1), size(H.V,1)) / 16)
            if ~isempty(nrm)                
                H = hmatrix_from_adaptive(@(i, j) H.U(i, :) * H.V(j, :)', ...
                        size(H.U,1), size(H.V,1), maxrank, @(i,j) 1, nrm);
            else
                H = hmatrix_from_adaptive(@(i, j) H.U(i, :) * H.V(j, :)', ...
                        size(H.U,1), size(H.V,1), maxrank);
            end
        end
    end
else
    % Recursive call on the children
    H.A11 = hmatrix_repartition(H.A11, maxrank, nrm);
    H.A12 = hmatrix_repartition(H.A12, maxrank, nrm);
    H.A21 = hmatrix_repartition(H.A21, maxrank, nrm);
    H.A22 = hmatrix_repartition(H.A22, maxrank, nrm);
    
    % Check if the blocks here are all full --- if they are, then we
    % can merge them and avoid having a deep tree for no gain.
    if ~isempty(H.A11.F) && ~isempty(H.A12.F) && ...
            ~isempty(H.A21.F) && ~isempty(H.A22.F)
        H.F = [ H.A11.F , H.A12.F ; H.A21.F , H.A22.F ];
        H.A11 = []; H.A12 = []; H.A21 = []; H.A22 = [];
    elseif H.A11.admissible && H.A12.admissible && H.A21.admissible && H.A22.admissible % if we have 4 children of low-rank let's try to merge them
        % first we compute the low-rank representation of the 4 blocks together
        [U, V] = merge_low_rank(H.A11, H.A12, H.A21, H.A22);
        
        if size(U, 2) < min([size(H) / 2, maxrank]) % if the rank is not too high...
            mem = sum(H.sz) * size(U, 2);
            mem11 = sum(H.A11.sz) * size(H.A11.U, 2);
            mem12 = sum(H.A12.sz) * size(H.A12.U, 2);
            mem21 = sum(H.A21.sz) * size(H.A21.U, 2);
            mem22 = sum(H.A22.sz) * size(H.A22.U, 2);
            
            if mem <= 1.75 * (mem11 + mem12 + mem21 + mem22) % ...and we do not lose too much in memory, then we merge them
                H.A11 = []; H.A12 = []; H.A21 = []; H.A22 = [];
                H.U = U;
                H.V = V;
                H.admissible = true;
            end
        end
    end
end

end

function [U, V] = merge_low_rank(A11, A12, A21, A22, nrm)
[U, V] = compress_factors(blkdiag([ A11.U, A12.U ], [ A21.U, A22.U ]), ...
    [ blkdiag(A11.V, A12.V), blkdiag(A21.V, A22.V) ]);
end

