function HU = hmatrix_chol(H)
%HMATRIX_CHOL Cholesky factorization for hm objects

HU = H;

if ~isempty(H.F)
    HU.F = chol(H.F);
else
    mp = size(H.A11,2);
    n = H.sz(2);
    HU.U21 = zeros(n-mp,0);
    HU.V21 = zeros(mp,0);
    HU.A11 = hmatrix_chol(H.A11);
    HU.U12 = HU.A11'\H.U12;
    HU.A22 = hmatrix_chol(hmatrix_rank_update(H.A22,-HU.V12*(HU.U12'*HU.U12),HU.V12));
end



end
