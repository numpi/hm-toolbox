function k = hss_rank(A)
if A.leafnode == 1
    k = 0;
else
    k = max([rank(A.B12), rank(A.B21), hss_rank(A.A11), hss_rank(A.A22)]);
end
end
