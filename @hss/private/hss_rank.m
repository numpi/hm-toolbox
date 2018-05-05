function k = hss_rank(A)
	if A.leafnode == 1
		k = 0;
	else
		k = max([min(size(A.B12)), min(size(A.B21)), hss_rank(A.A11), hss_rank(A.A22)]);
	end
end
