function H = hmatrix_mtimes(H1, H2)
%HMATRIX_MTIMES Multiply two H matrices. 

H = H1;

if ~isempty(H1.F)
	H.F = H1.F * H2.F;
else		
	H.A11 = hmatrix_rank_update(...
		hmatrix_mtimes(H1.A11, H2.A11), ...
		H1.U12 * (H1.V12' * H2.U21), H2.V21, true);
	
	H.A22 = hmatrix_rank_update(...
		hmatrix_mtimes(H1.A22, H2.A22), ...
		H1.U21 * (H1.V21' * H2.U12), H2.V12, true);
	
	[H.U12, H.V12] = compress_factors(...
		[ hmatrix_mtimes_dense(H1.A11, H2.U12), H1.U12 ], ...
		[ H2.V12, hmatrix_mtimes_dense(H2.A22', H1.V12) ]);
	
	[H.U21, H.V21] = compress_factors(...
		[ hmatrix_mtimes_dense(H1.A22, H2.U21), H1.U21 ], ...
		[ H2.V21, hmatrix_mtimes_dense(H2.A11', H1.V21) ]);
end

end

