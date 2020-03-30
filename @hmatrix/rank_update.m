function H = rank_update(H, U, V)
%RANK_UPDATE 
nrm = normest_Afun(@(x) H  * x + U * (V' * x), @(x) H' * x + V * (U' * x), size(H, 2), 1e-3);
H = hmatrix_rank_update(H, U, V, nrm);

end

