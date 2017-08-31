function k = hss_rank(A)
	if A.leafnode == 1
		k = 0;
	else
		k = max([size(A.Bu), size(A.Bl), hss_rank(A.hssl), hss_rank(A.hssr)]);
	end
end
