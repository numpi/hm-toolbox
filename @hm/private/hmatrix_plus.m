function H = hmatrix_plus(H1, H2)
%HMATRIX_PLUS Sum two hm objects

H = H1;

if ~isempty(H1.F)
	H.F = H1.F + H2.F;
else
	H.A11 = hmatrix_plus(H1.A11, H2.A11);
	H.A22 = hmatrix_plus(H1.A22, H2.A22);

	H.U12 = [ H1.U12, H2.U12 ];
	H.V12 = [ H1.V12, H2.V12 ];
	
	H.U21 = [ H1.U21, H2.U21 ];
	H.V21 = [ H1.V21, H2.V21 ];
end

end

