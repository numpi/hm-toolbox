function H = hmatrix_rank_update(H, U, V)
%HMATRIX_RANK_UPDATE Perform a low rank update to H. 

if ~isempty(H.F)
	H.F = H.F + U*V';
else
	mp = size(H.A11, 1);
	H.A11 = hmatrix_rank_update(H.A11, U(1:mp,:), V(1:mp,:));
	H.A22 = hmatrix_rank_update(H.A22, U(mp+1:end,:), V(mp+1:end,:));
		
	H.U21 = [ H.U21, U(mp+1:end,:) ];
	H.V21 = [ H.V21, V(1:mp,:) ];
	H.U12 = [ H.U12, U(1:mp,:) ];
	H.V12 = [ H.V12, V(mp+1:end,:) ];

end

