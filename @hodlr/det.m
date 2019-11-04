function d = det(H)
%DET Determinant of a HODLR matrix.
%
% D = DET(H) returns the determinant of H obtained by an LU factorization.

if size(H, 1) ~= size(H, 1)
    error('det(H): Matrix H is not square');
end

[L, U] = lu(H);

% The LU of the diagonal blocks is done using a pivoting strategy -- so we
% cannot assume that det(L) == 1
d = det_rec(U) * det_rec(L);

end

function d = det_rec(U)
if is_leafnode(U)
    if size(U, 1) ~= size(U, 2)
        error('det(H): Diagonal blocks need to be square');
    end
    
    d = det(U.F);
else
    d = det(U.A11) * det(U.A22);
end
end

