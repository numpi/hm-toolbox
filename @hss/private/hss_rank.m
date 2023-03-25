function k = hss_rank(A)
	k = hss_rank_rec(A, zeros(size(A.B12, 1), 0), zeros(size(A.B21, 1), 0));
end
function k = hss_rank_rec(A, M1, M2)
% Compute the rank of all the HSS block rows at each level
	if A.leafnode == 1
		k = 0;
	else
		k1 = rank([A.B12, M1]);
		k2 = rank([A.B21, M2]);
		if A.A11.leafnode == 0
			k3 = hss_rank_rec(A.A11, A.A11.Rl * [A.B12, M1], A.A11.Rr * [A.B12, M1]);
		else
			k3 = 0;
		end
		if A.A22.leafnode == 0	
			k4 = hss_rank_rec(A.A22, A.A22.Rl * [A.B21, M2], A.A22.Rr * [A.B21, M2]);
		else
			k4 = 0; 		
		end	
		k = max([k1, k2, k3, k4]);
	end	
end
