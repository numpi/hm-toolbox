function H = rank_update(H, U, V)
%RANK_UPDATE

% For the 2-norm we can take adavantage of the fast matrix-vector product
% for the low-rank part to estimate the norm beforehand. 
if hmatrixoption('norm') == 2
    nrm = normest_Afun(@(x) H  * x + U * (V' * x), ...
        @(x) H' * x + V * (U' * x), size(H, 2), 1e-3);
    H = hmatrix_rank_update(H, U, V, nrm);
else
    % Precomputing the Frobenius norm of H + UV' would help in this
    % process, as we could directly truncate relative to the norm of the
    % result.
    H = hmatrix_rank_update(H, U, V);
    H = hmatrix_compress(H);
end

end

