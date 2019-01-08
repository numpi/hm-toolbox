function IH = inv(H)
%INV computes the inverse of a HODLR matrix

if size(H, 1) ~= size(H, 2)
    error('Matrix is not square');
end

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
