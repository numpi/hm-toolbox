function k = hss_hodlr_rank(A)
% maximum rank of the off-diagonal blocks
	if A.leafnode == 1
		k = 0;
	else
		k = max([rank(A.B12), rank(A.B21), hss_hodlr_rank(A.A11), hss_hodlr_rank(A.A22)]);
	end
end
