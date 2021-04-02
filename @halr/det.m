function d = det(H)
%DET Determinant of a HALR matrix. 
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
	if U.sz(1) ~= U.sz(2)
		error('DET_REC:: non square diagonal block')
	end
	if U.admissible
		d = det(U.U * U.V');
        else
             	d = det(U.F);
        end
    else
        d = det(U.A11) * det(U.A22);
    end
end

