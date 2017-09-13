function nrm = norm_frobenius(A)
	A = hss_proper(A);
	nrm = sqrt(norm_frobenius_ric(A));
end
function nrm = norm_frobenius_ric(A)
	if A.leafnode == 1
		nrm = norm(A.D,'fro')^2;
	else
		nrm = sum(sum(A.Bu.^2) + sum(sum(A.Bl.^2)) + norm_frobenius_ric(A.hssl) + norm_frobenius_ric(A.hssr); 
	end
end
