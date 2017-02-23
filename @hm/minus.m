function H = minus(H1,H2)
%MINUS difference of matrices
if ~isa(H1,'hm') || ~isa(H2,'hm')
	error('Unsupported operation');
else
	H = H1;
	if ~isempty(H.F)
		H.F = H1.F - H2.F;
	else
		H.A11 = H1.A11 - H2.A11;
		H.A22 = H1.A22 - H2.A22;

		[H.U12, H.V12] = compress_factors([ H1.U12, -H2.U12 ], ...
									 [ H1.V12, H2.V12 ]);
		[H.U21, H.V21] = compress_factors([ H1.U21, -H2.U21 ], ...
									 [ H1.V21, H2.V21 ]);
	end
end
