function [HL, HU] = hmatrix_lu(H)
%HMATRIX_LU LU factorization for hm objects

HL = H;
HU = H;

if ~isempty(H.F)
	[HL.F, HU.F] = lu(H.F);
else
	mp = H.A11.sz(2);
	n = H.sz(2);
	HL.U12 = zeros(mp,0);
	HL.V12 = zeros(n-mp,0);
	HU.U21 = zeros(n-mp,0);
	HU.V21 = zeros(mp,0);
	[HL.A11, HU.A11] = hmatrix_lu(H.A11);
	HU.U12 = solve_lower_triangular(HL.A11, H.U12);
	% HL.V21 = (H.V21' / HU.A11)';
	HL.V21 = solve_upper_triangular2(HU.A11, H.V21')';
	[HL.A22, HU.A22] = hmatrix_lu(hmatrix_rank_update(H.A22,-HL.U21*(HL.V21'*HU.U12),HU.V12));
end


end

