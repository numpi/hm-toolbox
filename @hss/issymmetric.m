function y = issymmetric(A)
if A.leafnode == 1
	if ~issymmetric(A.D)
		y = false;
	else
		y = isequal(A.U, A.V);
	end
else
	if ~issymmetric(A.A11) || ~issymmetric(A.A22)
		y = false;
	else
		y = isequal(A.B12, A.B21') && isequal(A.Rl, A.Wl) && isequal(A.Rr, A.Wr);
	end
end
