function tH = tril(H)
%TRIL computes the lower triangular part of a HODLR matrix

if ~check_cluster_equality(H)
    error('tril is supported only for square matrices with square diagonal blocks');
end
tH = H;
if is_leafnode(H)
    tH.F = tril(H.F);
else
    tH.U12 = zeros(size(H.U12, 1), 0);
    tH.V12 = zeros(size(H.V12, 1), 0);
    tH.A11 = tril(H.A11);
    tH.A22 = tril(H.A22);
end
