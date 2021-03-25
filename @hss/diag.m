function d = diag(H)
%DIAG Obtain the diagonal of an HSS-matrix.

if H.leafnode
    d = diag(H.D);
else
    d = [ diag(H.A11) ; diag(H.A22) ];
end

end

