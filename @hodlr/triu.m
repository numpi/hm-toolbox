function tH = triu(H)
%TRIU computes the lower triangular part of a HODLR matrix

if ~check_cluster_equality(H)
    error('triu is supported only for square matrices with square diagonal blocks');
end
tH = H;
if is_leafnode(H)
    tH.F = triu(H.F);
else
    tH.U21 = zeros(size(H.U21, 1), 0);
    tH.V21 = zeros(size(H.V21, 1), 0);
    tH.A11 = triu(H.A11);
    tH.A22 = triu(H.A22);
end
