function H = cluster(H)
%CLUSTER Extract the tree representation of the cluster of the @hmatrix.

H.U = [];
H.V = [];
H.F = [];

if ~is_leafnode(H)
    H.A11 = cluster(H.A11);
    H.A12 = cluster(H.A12);
    H.A21 = cluster(H.A21);
    H.A22 = cluster(H.A22);
end

end

