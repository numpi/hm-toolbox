function IH = inv(H)
%INV computes the inverse of a HODLR matrix
IH = H;
if ~isempty(H.F)
    IH.F = inv(H.F);
elseif isempty(H.U12)
    IH = inv_lower_triangular(H);
elseif isempty(H.U21)
    IH = inv_upper_triangular(H);
else
    [HL, HU] = hmatrix_lu(H);
    IH = HU \ inv_lower_triangular(HL);
end
