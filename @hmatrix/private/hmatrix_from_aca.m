function H = hmatrix_from_aca(H, Afun, m, n)

if isempty(H)
	H = hmatrix_build_default_tree(m, n, hmatrixoption('block-size'));
end

H = hmatrix_from_aca_rec(H, Afun);

end

function H = hmatrix_from_aca_rec(H, Afun)
	[m, n] = size(H);
	if is_leafnode(H)
		if H.admissible
			[H.U, H.V] = aca(Afun, m, n, hmatrixoption('threshold'));
		else
			H.F = Afun((1:m).', (1:n));
		end
	else
		[m1, n1] = size(H.A11);
		H.A11 = hmatrix_from_aca_rec(H.A11, Afun);
		H.A12 = hmatrix_from_aca_rec(H.A12, @(i, j) Afun(i, j + n1));
		H.A21 = hmatrix_from_aca_rec(H.A21, @(i, j) Afun(i + m1, j));
		H.A22 = hmatrix_from_aca_rec(H.A22, @(i, j) Afun(i + m1, j + n1));
	end
end
